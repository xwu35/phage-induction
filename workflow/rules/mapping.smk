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

    rule genome_mapping_bowtie2:
        """map reads back to the genome using bowtie2"""
        input:
            genome_idx=multiext(os.path.join(dir["output"]["mapping"], "bowtie2_index", "{genome}"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}.sam"))
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
    rule genome_mapping_minimap2:
        """map reads back to genome using minimap2"""
        input:
            genome=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][1],
            R1=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}.sam"))
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
        sam=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}.sam")
    output:
        sorted=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}_sorted.bam"),
        index=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}_sorted.bam.bai")
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
    """generate stats"""
    input:
        sorted=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}_sorted.bam")
    output:
        depth=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}_depth.txt"),
        stats=os.path.join(dir["output"]["mapping"], "mapping_output", "{genome}", "{sample}_stats.txt")
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