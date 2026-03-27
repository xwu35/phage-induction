#------------ SET UP THE DIRECTORIES
dir = dict()
dir["output"] = dict()

# WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"] = os.path.join(workflow.basedir, "..", "db")

# OUTPUT DIRs
dir["output"]["base"] = RESULTS_DIR
dir["output"]["qc"] = os.path.join(dir["output"]["base"], "qc")
dir["output"]["trimmomatic"] = os.path.join(dir["output"]["qc"], "trimmomatic")
dir["output"]["fastqc"] = os.path.join(dir["output"]["qc"], "fastqc")
dir["output"]["reads_statistics"] = os.path.join(dir["output"]["qc"], "reads_statistics")
dir["output"]["mapping"] = os.path.join(dir["output"]["base"], "mapping")
dir["output"]["assembly"] = os.path.join(dir["output"]["base"], "assembly")
dir["output"]["alignment"] = os.path.join(dir["output"]["base"], "alignment")
dir["output"]["annotation"] = os.path.join(dir["output"]["base"], "annotation")

#------------ SET UP THE OUTPUT
pharokka=[]
for sample in SAMPLE:
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
mapped_within_prophage=[] # for counting primarily mapped reads within prophage region
for sample,genome in genome_name_dict.items():
    blastdb_index[sample]=[os.path.join(dir["output"]["alignment"], "blastDB", genome) + ext for ext in [".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"]]
    blastdb_prefix[sample]=os.path.join(dir["output"]["alignment"], "blastDB", genome)
    alignment.append(os.path.join(dir["output"]["alignment"], "blastn_filtered_within_prophage_region", genome, sample + "_blastn_filtered_within_prophage_region.tsv"))
    mapping.append(os.path.join(dir["output"]["mapping"], "mapping_stats", genome, sample + "_depth.txt"))
    mapping.append(os.path.join(dir["output"]["mapping"], "mapping_stats", genome, sample + "_stats.txt"))
    mapping.append(os.path.join(dir["output"]["mapping"], "bedgraph", genome, sample + ".bedgraph"))
    mapping.append(os.path.join(dir["output"]["mapping"], "read_counts", "prophage_region_counts", genome, sample + "_prophage_region_counts.tsv"))
    bowtie2_index[sample]=[os.path.join(dir["output"]["mapping"], "bowtie2_index", genome) + ext for ext in extensions]
    bowtie2_prefix[sample]=os.path.join(dir["output"]["mapping"], "bowtie2_index", genome)
    sam_files[sample]=os.path.join(dir["output"]["mapping"], "bam_files", genome, sample + ".sam")
    # for listing files as input
    mapped_bam_files.append(os.path.join(dir["output"]["mapping"], "bam_files", genome, sample + "_mapped_sorted.bam"))
    mapped_within_prophage.append(os.path.join(dir["output"]["mapping"], "bam_files", genome, sample + "_reads_within_prophage_region.bam"))

# raw reads fastqc
fastqc_input = [
    os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R2_stats.tsv"),
    os.path.join(dir["output"]["fastqc"], "multiqc_raw_reads", "multiqc_report.html")
]

# trimming
trimming_input = [
    os.path.join(dir["output"]["fastqc"], "multiqc_after_trimmomatic", "multiqc_report.html"),
    os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
    os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R2_stats.tsv")
    ] 

# mapping
mapping_input = mapping + [
    os.path.join(dir["output"]["mapping"], "read_counts", "overall", "all_mapped_read_counts.tsv"),
    os.path.join(dir["output"]["mapping"], "read_counts", "overall", "number_of_mapped_and_unmapped_reads.txt"),
    os.path.join(dir["output"]["mapping"], "read_counts", "overall", "reads_composition_barplot.svg"),
    os.path.join(dir["output"]["mapping"], "read_counts", "overall", "mapped_read_counts_within_prophage_regions.txt")
]

# assembly
assembly_input = [
    os.path.join(dir["output"]["assembly"], "contig_evaluation", "combined_assembly_stats.txt"),
    os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")
    ]

# alignment
alignment_input = alignment

# annotation 
annotation_input = pharokka

# all steps excep fastqc
all_input = trimming_input + mapping_input + assembly_input + alignment_input + annotation_input

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