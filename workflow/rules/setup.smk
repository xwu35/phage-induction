#------------ SET UP THE DIRECTORIES
dir = dict()
dir["output"] = dict()

# WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"] = os.path.join(workflow.basedir, "..", "db")

# OUTPUT DIRs
dir["output"]["base"] = RESULTS_DIR
dir["output"]["trimmomatic"] = os.path.join(dir["output"]["base"], "trimmomatic")
dir["output"]["fastqc"] = os.path.join(dir["output"]["base"], "fastqc")
dir["output"]["reads_statistics"] = os.path.join(dir["output"]["base"], "reads_statistics")
dir["output"]["intermediate"] = os.path.join(dir["output"]["base"], "intermediate")
dir["output"]["mapping"] = os.path.join(dir["output"]["base"], "mapping")
dir["output"]["assembly"] = os.path.join(dir["output"]["base"], "assembly")
dir["output"]["alignment"] = os.path.join(dir["output"]["base"], "alignment")
dir["output"]["annotation"] = os.path.join(dir["output"]["base"], "annotation")
dir["output"]["viral_identification"] = os.path.join(dir["output"]["base"], "viral_identification")

#------------ SET UP THE OUTPUT
viral_identification=[]
pharokka=[]
for sample in SAMPLE:
    viral_identification.append(os.path.join(dir["output"]["viral_identification"], "blastn", sample, sample + "_blastn.out"))
    pharokka.append(os.path.join(dir["output"]["annotation"], "pharokka", "all_gff", sample + ".gff"))

# extract the first item of a list (genome name) for each key/sample in a dictionary
genome_name_dict = {key: value[0] for key, value in SAMPLE_GENOME_MAP.items()}

blastdb_index={}
blastdb_prefix={}
alignment=[]
mapping=[]
bowtie2_index={}
bowtie2_prefix={}
extensions = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
sam_files={}
mapped_bam_files=[] # for counting all primary reads
mapped_prophage_files=[] # for counting primarily mapped reads to prophage region
mapped_within_prophage=[] # for counting primarily mapped reads within prophage region
for sample,genome in genome_name_dict.items():
    blastdb_index[sample]=[os.path.join(dir["output"]["alignment"], "blastDB", genome) + ext for ext in [".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"]]
    blastdb_prefix[sample]=os.path.join(dir["output"]["alignment"], "blastDB", genome)
    alignment.append(os.path.join(dir["output"]["alignment"], "blastn_filtered_within_prophage_region", genome, sample + "_blastn_filtered_within_prophage_region.tsv"))
    alignment.append(os.path.join(dir["output"]["alignment"], "blastn", genome, sample + "_blastn_anicalc_selected.out"))
    mapping.append(os.path.join(dir["output"]["mapping"], "mapping_output", genome, sample + "_depth.txt"))
    mapping.append(os.path.join(dir["output"]["mapping"], "mapping_output", genome, sample + "_stats.txt"))
    mapping.append(os.path.join(dir["output"]["mapping"], "bedgraph", genome, sample + ".bedgraph"))
    mapping.append(os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", genome, sample + "_counts.txt"))
    mapping.append(os.path.join(dir["output"]["mapping"], "region_counts", genome, sample + "_region_counts.tsv"))
    bowtie2_index[sample]=[os.path.join(dir["output"]["mapping"], "bowtie2_index", genome) + ext for ext in extensions]
    bowtie2_prefix[sample]=os.path.join(dir["output"]["mapping"], "bowtie2_index", genome)
    sam_files[sample]=os.path.join(dir["output"]["mapping"], "mapping_output", genome, sample + ".sam")
    # for listing files as input
    mapped_bam_files.append(os.path.join(dir["output"]["mapping"], "filtered_bam", genome, sample + "_mapped_sorted.bam"))
    mapped_prophage_files.append(os.path.join(dir["output"]["mapping"], "filtered_bam", genome, sample + "_mapped_sorted_prophage.bam"))
    mapped_within_prophage.append(os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", genome, sample + "_contained_reads.bam"))

# checkv scripts
scripts_input = [
    os.path.join(dir["scripts"], "anicalc.py"),
    os.path.join(dir["scripts"], "aniclust.py")
]

# raw reads fastqc
fastqc_input = [
    os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R2_stats.tsv"),
    os.path.join(dir["output"]["fastqc"], "multiqc_raw_reads", "multiqc_report.html")
]

# trimming
trimming_input = scripts_input + [
    os.path.join(dir["output"]["fastqc"], "multiqc_after_trimmomatic", "multiqc_report.html"),
    os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R2_stats.tsv")
    ] 

# mapping
mapping_input = mapping + [
    os.path.join(dir["output"]["mapping"], "read_counts", "mapped_reads_counts.tsv"),
    os.path.join(dir["output"]["mapping"], "read_counts", "all_prophage_regions_read_counts.txt"),
    os.path.join(dir["output"]["mapping"], "read_counts", "within_prophage_regions_read_counts.txt")
]

# assembly
assembly_input = [
    os.path.join(dir["output"]["assembly"], "contig_evaluation", "combined_assembly_stats.txt"),
    os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")
    ]

# identification
identification_input = viral_identification

# alignment
alignment_input = scripts_input + alignment

# annotation 
annotation_input = pharokka

#---------- DOWNLOAD SCRIPTS
localrules:
    download_checkV_scripts

rule download_checkV_scripts:
    """
    download scripts for calculating ANI and clustering
    """
    output:
        anicalc=os.path.join(dir["scripts"], "anicalc.py"),
        aniclust=os.path.join(dir["scripts"], "aniclust.py")
    params:
        dir=dir["scripts"]
    shell:
        """
        # download the two scripts from checkV which cannot be installed via conda
        # only download them if the files don't exsit
        # script for calculating ANI and AF
        if [[ ! -e {output.anicalc} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/anicalc.py && chmod +x {output.anicalc}
        else
            # in case the files were downloaded, but not executable 
            chmod +x {output.anicalc}
        fi

        # script for clustering contigs into vOTUs
        if [[ ! -e {output.aniclust} ]]; then 
            wget -P {params.dir} https://bitbucket.org/berkeleylab/checkv/raw/3f185b5841e8c109848cd0b001df7117fe795c50/scripts/aniclust.py && chmod +x {output.aniclust}
        else
            # in case the files were downloaded, but not executable
            chmod +x {output.aniclust}
        fi
        """

#---------- DOWNLOAD DATABASE
rule download_pharokka_db:
    output:
        os.path.join(dir["db"], "pharokka_db", ".done")
    params:
        outdir=os.path.join(dir["db"], "pharokka_db")
    threads:
        config["resources"]["small_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "pharokka.yml")
    shell:
        """
        install_databases.py -o {params.outdir} &&
        touch {output}
        """