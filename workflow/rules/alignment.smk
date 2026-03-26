rule create_blastdb:
    """create blast database for bacterial genome"""
    input:
        lambda wildcards: GENOME_SEQUENCE_MAP[wildcards.genome]
    output:
        blastdb=multiext(os.path.join(dir["output"]["alignment"], "blastDB", "{genome}"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
    params: 
        prefix=os.path.join(dir["output"]["alignment"], "blastDB", "{genome}"),
        dbtype=config["blastn"]["dbtype"]
    threads: 
        config["resources"]["small_cpu"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        makeblastdb -in {input} -dbtype {params.dbtype} -out {params.prefix} 
        """

rule blastn:
    """blast against the bacterial genome"""
    input:
        blastdb=multiext(os.path.join(dir["output"]["alignment"], "blastDB", "{genome}"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        blastn=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn.out")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    params:
        db_path=os.path.join(dir["output"]["alignment"], "blastDB", "{genome}"),
        evalue=config["blastn"]["evalue"],
        max_target_seqs=config["blastn"]["max_target_seqs"],
        outfmt=config["blastn"]["outfmt"]
    conda:
        os.path.join(dir["env"], "alignment.yml")
    shell:
        """
        blastn -db {params.db_path} \
            -query {input.contigs} \
            -out {output.blastn} \
            -outfmt {params.outfmt} \
            -max_target_seqs {params.max_target_seqs} \
            -num_threads {threads} \
            -evalue {params.evalue} 
        """

rule filter_blastn_alignment:
    """
    keep contigs with >=95% sequences aligned to genome with >=98% identity
    """
    input:
        blastn=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn.out")
    output:
        filtered=os.path.join(dir["output"]["alignment"], "blastn_filtered", "{genome}", "{sample}_blastn_filtered.tsv")
    params:
        script=os.path.join(dir["scripts"], "filter_blastn_alignments.py")
    shell:
        """
        # if blastn output is not empty, then filter
        if [[ -s {input.blastn} ]]; then
            {params.script} -i {input.blastn} -o {output.filtered}
        else
            echo "No contigs aligned to the bacterial genome. BLASTn result is empty."
            echo "chr\tstart\tend" > {output.filtered} # add header line for visualization app
        fi
        """

rule get_blastn_alignment_within_prophage_region:
    """
    keep contigs aligned within prophage region
    """
    input:
        filtered=os.path.join(dir["output"]["alignment"], "blastn_filtered", "{genome}", "{sample}_blastn_filtered.tsv"),
        bed_file=lambda wildcards: SAMPLE_GENOME_MAP[wildcards.sample][2]
    output:
        os.path.join(dir["output"]["alignment"], "blastn_filtered_within_prophage_region", "{genome}", "{sample}_blastn_filtered_within_prophage_region.tsv")
    params:
        script=os.path.join(dir["scripts"], "get_alignments_within_prophage_region.py")
    shell:
        """
        # if filtered blastn output is not empty, then filter
        if [[ -s {input.filtered} ]]; then
            {params.script} -b {input.filtered} -p {input.bed_file} -o {output}
        else
            echo "The filtered BLASTn file is empty."
            echo "chr\tstart\tend\tcontig\tstart_p\tend_p\tregion" >{output} # add header line for visualization app
        fi
        """

rule blastn_ani_calculation:
    """
    caculate ANI from blastn
    """
    input:
        blastn=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn.out")
    output:
        ani=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn_anicalc.out")
    params:
        anicalc=os.path.join(dir["scripts"], "anicalc.py")
    shell:
        """
        # if blastn output is not empty, then calculate ANI
        if [[ -s {input.blastn} ]]; then
            {params.anicalc} -i {input.blastn} -o {output.ani}
        else
            echo "No contigs aligned to the bacterial genome. BLASTn result is empty."
            touch {output.ani}
        fi
        """

rule blastn_contig_selection:
    """
    select the contigs aligned to the bacterial genome with ANI >= 95% and AF >=85%
    """
    input:
        ani=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn_anicalc.out")
    output:
        selected=os.path.join(dir["output"]["alignment"], "blastn", "{genome}", "{sample}_blastn_anicalc_selected.out")
    shell:
        """
        # if ani is not empty, then select contigs based on ANI and AF
        if [[ -s {input.ani} ]]; then
            awk '{{OFS="\t"}}{{FS="\t"}}{{ if (NR==1) {{print $1, $2, $4, $5}} else if ($4>=95 && $5>=85) {{print $1, $2, $4, $5}}}}' {input.ani} > {output.selected}
        else
            echo "No contigs aligned to the bacterial genome. ANI File is empty"
            touch {output.selected}
        fi
        """
