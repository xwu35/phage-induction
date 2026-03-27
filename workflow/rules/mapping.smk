if config["mapper"]=="bowtie2":
    rule bowtie2_build_genome:
        """build bowtie2 index for the genome"""
        input:
            lambda wildcards: GENOME_SEQUENCE_MAP[wildcards.genome]
        output:
            multiext(os.path.join(dir["output"]["mapping"], "bowtie2_index", "{genome}"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=os.path.join(dir["output"]["mapping"], "bowtie2_index", "{genome}")
        threads: 
            config["resources"]["small_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input} \
                {params.prefix} 
            """

    rule bowtie2_mapping:
        """map reads back to the genome using bowtie2"""
        input:
            genome_idx=multiext(os.path.join(dir["output"]["mapping"], "bowtie2_index", "{genome}"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=lambda wildcards: bowtie2_prefix[wildcards.sample]
        log:
            os.path.join(os.path.join(dir["output"]["mapping"], "logs", "{genome}", "{sample}.bowtie2.log"))
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} -1 {input.R1} -2 {input.R2} \
                -S {output.sam} 2> {log} 
            """

elif config["mapper"]=="minimap2":
    rule minimap2_mapping:
        """map reads back to genome using minimap2"""
        input:
            genome=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][1],
            R1=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"]
        log:
            os.path.join(dir["output"]["mapping"], "logs", "{genome}", "{sample}.minimap2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no \
                {input.genome} \
                {input.R1} \
                {input.R2} > {output.sam} 2>{log}
            """

rule sort_and_index_bam:
    """sort and index bam file"""
    input:
        sam=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.sam")
    output:
        sorted=temp(os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.s.bam")),
        index=temp(os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.s.bam.bai"))
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools view -bS {input.sam} | \
        samtools sort -@ {threads} -T /tmp -o {output.sorted}
        samtools index {output.sorted} -o {output.index}
        """

rule generate_stats:
    input:
        sorted=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.s.bam"),
        index=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.s.bam.bai")
    output:
        depth=os.path.join(dir["output"]["mapping"], "mapping_stats", "{genome}", "{sample}_depth.txt"),
        stats=os.path.join(dir["output"]["mapping"], "mapping_stats", "{genome}", "{sample}_stats.txt")
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools depth -aa {input.sorted} > {output.depth}
        samtools stats {input.sorted} > {output.stats}
        """

rule keep_mapped_reads:
    """keep mapped reads"""
    input:
        sam=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}.sam")
    output:
        mapped=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_mapped_sorted.bam"),
        index=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_mapped_sorted.bam.bai")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # keep only mapped reads and sort the bam file by reference coordinates 
        samtools view -bS -F 4 {input.sam} | samtools sort - -o {output.mapped} -@ {threads}

        # index bam file
        samtools index {output.mapped} -o {output.index}
        """

rule count_all_primarily_mapped_raeds:
    input:
        mapped_bam_files
    output:
        os.path.join(dir["output"]["mapping"], "read_counts", "overall", "all_mapped_read_counts.tsv")
    params:
        dir=os.path.join(dir["output"]["mapping"], "bam_files"),
        script=os.path.join(dir["scripts"], "count_primary_alignments.sh")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        {params.script} -i {params.dir} -e _mapped_sorted.bam -o {output}
        """

rule plot_reads_distribution:
    """combine reads counts and calculate number of removed reads at each step"""
    input:
        trimmed=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
        mapped=os.path.join(dir["output"]["mapping"], "read_counts", "overall", "all_mapped_read_counts.tsv")
    output:
        table=os.path.join(dir["output"]["mapping"], "read_counts", "overall", "number_of_mapped_and_unmapped_reads.txt"),
        figure=os.path.join(dir["output"]["mapping"], "read_counts", "overall", "reads_composition_barplot.svg")
    params:
        script=os.path.join(dir["scripts"], "combine_read_counts_and_plot.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    script:
        "{params.script}"

rule generate_coverage_bedgraph:
    input:
        os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_mapped_sorted.bam")
    output:
        os.path.join(dir["output"]["mapping"], "bedgraph", "{genome}", "{sample}.bedgraph")
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """    
        bedtools genomecov -bga -ibam {input} > {output}
        """

rule keep_reads_within_prophage_region:
    input:
        bam=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_mapped_sorted.bam"),
        index=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_mapped_sorted.bam.bai"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_reads_within_prophage_region.bam")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "bedtools.yml")
    shell:
        """
        # -f 1: the read must overlap the region by at least 1.0 (100%) of its total length
        bedtools intersect -f 1 -wa -a {input.bam} -b {input.bed_file} > {output}
        """

rule count_primarily_mapped_raeds_within_prophage_region:
    input:
        mapped_within_prophage
    output:
        os.path.join(dir["output"]["mapping"], "read_counts", "overall", "mapped_read_counts_within_prophage_regions.txt")
    params:
        dir=os.path.join(dir["output"]["mapping"], "bam_files"),
        script=os.path.join(dir["scripts"], "count_primary_alignments.sh")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        {params.script} -i {params.dir} -e _reads_within_prophage_region.bam -o {output}
        """

rule count_mapped_reads_for_each_prophage_region:
    """count number of reads for each prophage region"""
    input:
        bam=os.path.join(dir["output"]["mapping"], "bam_files", "{genome}", "{sample}_reads_within_prophage_region.bam"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        temp(os.path.join(dir["output"]["mapping"], "read_counts", "prophage_region_counts", "{genome}", "{sample}_counts.txt"))
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
        region=os.path.join(dir["output"]["mapping"], "read_counts", "prophage_region_counts", "{genome}", "{sample}_counts.txt"),
        trimmed=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv")
    output:
        os.path.join(dir["output"]["mapping"], "read_counts", "prophage_region_counts", "{genome}", "{sample}_prophage_region_counts.tsv")
    params:
        script=os.path.join(dir["scripts"], "calculate_proportion_reads.py")
    shell:
        """    
        {params.script} -i {input.region} -t {input.trimmed} -o {output}
        """
        