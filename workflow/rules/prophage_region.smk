rule filter_out_unmapped_reads:
    """filter out unmapped reads"""
    input:
        lambda wildcards: sam_files[wildcards.sample]
    output:
        mapped=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # keep only mapped reads and sort the bam file by reference coordinates 
        # -F 4 flag to filter out unmapped reads
        samtools view -bS -F 4 {input} | samtools sort - -o {output.mapped} -@ {threads}
        """

rule count_all_primarily_mapped_raeds:
    input:
        mapped_bam_files
    output:
        os.path.join(dir["output"]["mapping"], "read_counts", "mapped_reads_counts.tsv")
    params:
        dir=os.path.join(dir["output"]["mapping"], "filtered_bam"),
        script=os.path.join(dir["scripts"], "count_primary_alignments.sh")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        {params.script} -i {params.dir} -e _mapped_sorted.bam -o {output}
        """

rule combine_all_mapped_reads_status:
    """combine reads counts and calculate number of removed reads at each step NEED UPDATE"""
    input:
        trimmed=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
        mapped=os.path.join(dir["output"]["mapping"], "read_counts", "mapped_reads_counts.tsv")
    output:
        table=os.path.join(dir["output"]["mapping"], "read_counts", "combined", "number_of_mapped_and_unmapped_reads.txt"),
        figure=os.path.join(dir["output"]["mapping"], "read_counts", "combined", "reads_composition_barplot.svg")
    params:
        script=os.path.join(dir["scripts"], "combine_read_counts_and_plot.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    script:
        "{params.script}"

rule host_separate_prophage_and_nonphage:
    """separate prophage and nonprophage regions"""
    input:
        mapped=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted.bam"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        prophage=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_prophage.bam"),
        prophage_index=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_prophage.bam.bai"),
        nonprophage=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_non-prophage.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # index the bam file
        samtools index {input.mapped} 

        # separate mapped prophage and non-prophage reads
        # {output.prophage} will contain alignments that overlap with regions in the bed file
        samtools view {input.mapped} \
            -b -h -o {output.prophage} \
            -U {output.nonprophage} \
            -L {input.bed_file} \
            -@ {threads}
        
        samtools index {output.prophage} -o {output.prophage_index}
        """

rule generate_coverage_file:
    input:
        os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_prophage.bam")
    output:
        os.path.join(dir["output"]["mapping"], "bedgraph", "{genome}", "{sample}.bedgraph")
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """    
        bedtools genomecov -bga -ibam {input} > {output}
        """

rule only_keep_alignment_within_region:
    input:
        bam=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_prophage.bam"),
        index=os.path.join(dir["output"]["mapping"], "filtered_bam", "{genome}", "{sample}_mapped_sorted_prophage.bam.bai"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", "{genome}", "{sample}_contained_reads.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """
        # -f 1: the read must overlap the region by at least 1.0 (100%) of its total length
        bedtools intersect -f 1 -wa -a {input.bam} -b {input.bed_file} > {output}
        """

rule count_reads_within_region:
    """count number of reads for each prophage region"""
    input:
        bam=os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", "{genome}", "{sample}_contained_reads.bam"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", "{genome}", "{sample}_counts.txt")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """    
        bedtools coverage -counts -a {input.bed_file} -b {input.bam} > {output}
        """

rule calculate_proportion_reads_per_region:
    """calculate the proportion of reads mapped to each region"""
    input:
        region=os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions", "{genome}", "{sample}_counts.txt"),
        trimmed=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
    output:
        os.path.join(dir["output"]["mapping"], "region_counts", "{genome}", "{sample}_region_counts.tsv")
    params:
        script=os.path.join(dir["scripts"], "calculate_proportion_reads.py")
    shell:
        """    
        {params.script} -i {input.region} -t {input.trimmed} -o {output}
        """
        
rule count_primary_alignments_at_prophage_region:
    input:
        prophage=mapped_prophage_files,
        within_regions=mapped_within_prophage
    output:
        all=os.path.join(dir["output"]["mapping"], "read_counts", "all_prophage_regions_read_counts.txt"),
        within=os.path.join(dir["output"]["mapping"], "read_counts", "within_prophage_regions_read_counts.txt")
    params:
        all_prophage_dir=os.path.join(dir["output"]["mapping"], "filtered_bam"),
        within_regions_dir=os.path.join(dir["output"]["mapping"], "alignments_within_prophage_regions"),
        script=os.path.join(dir["scripts"], "count_primary_alignments.sh")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        {params.script} -i {params.all_prophage_dir} -e _mapped_sorted_prophage.bam -o {output.all}
        {params.script} -i {params.within_regions_dir} -e _contained_reads.bam -o {output.within}
        """