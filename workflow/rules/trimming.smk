localrules:
    rename_raw_reads

rule rename_raw_reads:
    """rename raw reads to make sure the extension works for fastqc"""
    input:
        R1=lambda wildcards: R1_MAP[wildcards.sample],
        R2=lambda wildcards: R2_MAP[wildcards.sample]
    output:
        R1=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    shell:
        """
        ln -s {input.R1} {output.R1} 
        ln -s {input.R2} {output.R2} 
        """

rule fastqc_raw_reads:
    """run fastqc on raw reads"""
    input:
        R1=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["fastqc"], "raw_reads", "{sample}_R1_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["fastqc"], "raw_reads", "{sample}_R2_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1} -o {params.dir} -t {threads}
        fastqc {input.R2} -o {params.dir} -t {threads}
        """

rule multiqc_raw_reads:
    """aggregate fastqc results from raw reads"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["fastqc"], "raw_reads", "{sample}_R1_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["fastqc"], "raw_reads", "{sample}_R2_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["fastqc"], "multiqc_raw_reads", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["fastqc"], "raw_reads")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule get_raw_reads_status:
    """check raw reads counts"""
    input:
        R1=expand(os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R1.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R2.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_statistics"], "raw_reads", "R2_stats.tsv")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """

rule trimmomatic:
    """remove adapters and low quality reads"""
    input:
        R1=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["qc"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    output:
        R1P=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
        R1U=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.unpaired.fastq.gz"),
        R2P=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz"),
        R2U=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.unpaired.fastq.gz")
    params:
        seqType = config["trimmomatic"]["seqType"],
        phred = config["trimmomatic"]["phred"],
        adapter = config["adapter"],
        adapter_params = config["trimmomatic"]["adapter_params"],
        post_adapter_params = config["trimmomatic"]["post_adapter_params"],
    log:
        os.path.join(dir["output"]["trimmomatic"], "{sample}.trimmomatic.log"),
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "trimmomatic.yml")
    shell:
        """
        trimmomatic {params.seqType} {params.phred} -threads {threads} \
            {input.R1} {input.R2} \
            {output.R1P} {output.R1U} {output.R2P} {output.R2U} \
            ILLUMINACLIP:{params.adapter}:{params.adapter_params} \
            {params.post_adapter_params} 2>{log} # –baseout
        """

rule fastqc_trimmed_reads:
    """run fastqc on trimmed reads after trimmomatic to check if adapters were removed"""
    input:
        R1P=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
        R2P=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["fastqc"], "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["fastqc"], "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1P} -o {params.dir} -t {threads}
        fastqc {input.R2P} -o {params.dir} -t {threads}
        """

rule multiqc_trimmed_reads:
    """aggregate fastqc results from trimmed reads after trimmomatic"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["fastqc"], "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["fastqc"], "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["fastqc"], "multiqc_after_trimmomatic", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["fastqc"], "after_trimmomatic")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule get_reads_status_after_trimmomatic:
    """check reads counts after removing adapters and low quality reads"""
    input:
        R1=expand(os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"), sample=SAMPLE),
        R2=expand(os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R1_stats.tsv"),
        R2_stats=os.path.join(dir["output"]["reads_statistics"], "after_trimmomatic", "R2_stats.tsv")
    threads:
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        seqkit stats -j {threads} -To {output.R2_stats} {input.R2}
        """
