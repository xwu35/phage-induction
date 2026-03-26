rule spades_assembly:
    input:
        R1=os.path.join(dir["output"]["trimmomatic"], "{sample}_R1.trimmomatic.fastq.gz"),
        R2=os.path.join(dir["output"]["trimmomatic"], "{sample}_R2.trimmomatic.fastq.gz")
    output:
        contigs=os.path.join(dir["output"]["assembly"], "intermediate", "{sample}", "contigs.fasta"),
        renamed=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs.fasta")
    params:
        outdir=directory(os.path.join(dir["output"]["assembly"], "intermediate", "{sample}")),
        setting=config["spades"]["setting"]
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "spades.yml")
    shell:
        """
        # run SPADES
        spades.py {params.setting} \
            -1 {input.R1} \
            -2 {input.R2} \
            -o {params.outdir} \
            -t {threads} \
            --tmp-dir /tmp

        # rename contigs using the sample name
        sed 's/>/>{wildcards.sample}_/' {output.contigs} > {output.renamed}
        """

rule extract_1kb_contigs:
    """keep contigs >= 1kb"""
    input:
        renamed=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs.fasta")
    output:
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if contig file is not empty, then select contigs >= 1kb
        if [[ -s {input.renamed} ]]; then
            seqkit seq -m 1000 {input.renamed} > {output.contigs_1kb}
        else
            # if samples have no assembled contigs, then create an empty file to avoid error
            echo "No contigs were assembled from this sample"
            touch {output.contigs_1kb}
        fi
        """

rule assembly_status:
    """check the quality of contigs using Quast"""
    input:
        contigs_1kb=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        report=os.path.join(dir["output"]["assembly"], "contig_evaluation", "{sample}", "transposed_report.tsv")
    threads: 
        config["resources"]["small_cpu"]
    params:
        dir=lambda w, output: os.path.dirname(output.report),
        length=config["quast"]["length"]
    conda:
        os.path.join(dir["env"], "quast.yml")
    shell:
        """
        # if contigs_1kb is not empty, then evaluate contig status
        if [[ -s {input.contigs_1kb} ]]; then
            quast.py {input.contigs_1kb} -m {params.length} -t {threads} -o {params.dir}
        else
            echo "No contigs >= 1kb were obtained from this sample"
            touch {output.report}
        fi
        """
        
rule combine_assembly_status:
    """combine all quast assembly stats"""
    input:
        report=expand(os.path.join(dir["output"]["assembly"], "contig_evaluation", "{sample}", "transposed_report.tsv"), sample=SAMPLE)
    output:
        os.path.join(dir["output"]["assembly"], "contig_evaluation", "combined_assembly_stats.txt")
    params:
        dir=directory(os.path.join(dir["output"]["assembly"], "contig_evaluation")),
        script=os.path.join(dir["scripts"], "combine_quast_assembly_stats.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    shell:
        """
        Rscript {params.script} -d {params.dir} -o {output}
        """

rule combine_all_contigs_1kb:
    """
    combine contigs from all samples
    """
    input:
        expand(os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta"), sample=SAMPLE)
    output:
        os.path.join(dir["output"]["assembly"], "all_contigs_1kb.fasta")
    shell:
        """
        cat {input} > {output}
        """