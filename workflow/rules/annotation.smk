localrules:
    pool_gffs

rule pharokka_annotation:
    """functional annotation on the assembled contigs"""
    input:
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta"),
        done=os.path.join(dir["db"], "pharokka_db", ".done")
    output:
        blastn=os.path.join(dir["output"]["annotation"], "pharokka", "{sample}", "{sample}.gff")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    params:
        dir=os.path.join(dir["output"]["annotation"], "pharokka", "{sample}"),
        database=os.path.join(dir["db"], "pharokka_db"),
        meta_setting=config["pharokka"]["meta_setting"],
        regular_setting=config["pharokka"]["regular_setting"]
    conda:
        os.path.join(dir["env"], "pharokka.yml")
    shell:
        """
        # --meta mode doesn't work if the sample has only one contig
        count=$( grep '>' {input.contigs} | wc -l)

        # if number of contig > 1, then use --meta
        if [ $count -gt 1 ]; then
            pharokka.py \
                -i {input.contigs} \
                -o {params.dir} \
                -d {params.database} \
                -t {threads} -f -p {wildcards.sample} \
                {params.meta_setting}
        else
            pharokka.py \
                -i {input.contigs} \
                -o {params.dir} \
                -d {params.database} \
                -t {threads} -f -p {wildcards.sample} \
                {params.regular_setting}

        fi
        """

rule pool_gffs:
    """
    remove input sequences (deletes everything from the ##FASTA line to the end of the file),
    remove lines in the GFF file that start with ##,
    pool them into the same folder
    """
    input:
        os.path.join(dir["output"]["annotation"], "pharokka", "{sample}", "{sample}.gff")
    output:
        os.path.join(dir["output"]["annotation"], "pharokka", "all_gff", "{sample}.gff")
    shell:
        """
        cat {input} | sed '/^##FASTA$/,$d' | grep -v "^##" > {output}
        """