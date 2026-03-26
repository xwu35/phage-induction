rule genomad:
    """
    Viral contig identification using geNomad
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        genomad=directory(os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "genomad")),
        summary=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "genomad", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus_summary.tsv"),
        genomad_viruses=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "genomad", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus.fna")
    params:
        database=config["genomad"]["database_path"],
        sensitivity=config["genomad"]["sensitivity"], 
        other=config["genomad"]["other"]
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["big_mem"]
    shell:
        """
        {GENOMAD}

        if [ -s "{input.contigs}" ]; then
            genomad end-to-end \
                --cleanup {input.contigs} \
                {output.genomad} \
                {params.database} \
                -t {threads} \
                --sensitivity {params.sensitivity} \
                {params.other}
        else
            echo "File is empty"
            touch {output.summary}
            touch {output.genomad_viruses} # this one is necessary for extracting genomad identified provirus
        fi   
        """

rule virsorter2:
    """
    Viral contig identification using VirSorter2
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        virsorter2_outdir=directory(os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "virsorter2-pass1")),
        score=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "virsorter2-pass1", "final-viral-score.tsv")
    params:
        groups=config["virsorter2"]["groups"],
        length=config["virsorter2"]["length"],
        score=config["virsorter2"]["score"],
        other=config["virsorter2"]["other"],
        steps=config["virsorter2"]["steps"]
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["big_mem"]
    shell:
        """
        {VIRSORTER2}

        if [ -s "{input.contigs}" ]; then
            virsorter run \
                {params.other} \
                -i {input.contigs} \
                -w {output.virsorter2_outdir} \
                --include-groups {params.groups} \
                --min-length {params.length} \
                --min-score {params.score} \
                -j {threads} {params.steps}
        else
            echo "Contig file is empty"
            touch {output.score}
        fi   
        """

rule cenote_taker3:
    """
    Viral contig identification using Cenote-Taker3.
    """
    input: 
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        # cenote-taker3 doesn't output summary files if no viruses were found, so cannot track summary files
        done=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "cenote_taker3", ".done")
    params:
        # because snakemake generates the directory with sample name, cenote-taker3 will rename it (e.g. change VLP_Day0 to VLP_Day0_old_OJ5RL)
        # so here use a tmp directory and then move everything to destination directory
        tmp_dir=directory(os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "cenote_taker3_tmp")),
        dst_dir=directory(os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "cenote_taker3")),
        database=config["cenote_taker3"]["database"],
        # because the {} in -db {virion,rdrp,dnarep} makes snakemake thinking they are wildcards, 
        # therefore need to deactivate automatic wildcard expansion in params strings by specifying `lambda wildcards: `
        settings=lambda wildcards: config["cenote_taker3"]["settings"] 
    threads:
        config["resources"]["med_cpu"]
    resources:  
        mem_mb=config["resources"]["small_mem"]
    shell:
        """
        {CT3}

        if [ -s "{input.contigs}" ]; then
            cenotetaker3 -c {input.contigs} \
                -r cenote_taker3 \
                -t {threads} \
                -wd {params.tmp_dir} \
                --cenote-dbs {params.database} \
                {params.settings} && # && lets you do something based on whether the previous command completed successfully
            mv {params.tmp_dir}/cenote_taker3/* {params.dst_dir} &&
            rm -r {params.tmp_dir} &&
            touch {output.done}
        else
            echo "Contig file is empty"
            touch {output.done}
        fi   
        """

rule select_viral_contigs:
    """
    Contigs satisfying one of the following criteria were considered viral:
    1. geNomad: all identified contigs
    2. Cenote-taker3: all identified contigs
    3. VirSorter2: all identified contigs
    """
    input:
        genomad=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "genomad", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus_summary.tsv"),
        virsorter2=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "virsorter2-pass1", "final-viral-score.tsv"),
        done=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "cenote_taker3", ".done")
    output:
        viral=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_viral_contig_names.tsv"),
        provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_genomad_provirus_contig_names.tsv")
    params:
        script=os.path.join(dir["scripts"], "baldridge_select_viral_contigs.R"),
        cenote=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "cenote_taker3", "cenote_taker3_virus_summary.tsv")
    conda:
        os.path.join(dir["env"], "R.yml")
    shell:
        """
        # if Cenote-Taker3's summary file doesn't exist, create one before running the Rscript
        # Rscript deals with empty input file
        if [[ ! -f {params.cenote} ]]; then
            touch {params.cenote}
        fi 
        
        Rscript {params.script} \
                -g {input.genomad} \
                -c {params.cenote} \
                -v {input.virsorter2} \
                -o {output.viral} \
                -p {output.provirus}
        """

rule extract_identified_viral_contig:
    """
    extract_identified_viral_contig sequences. 
    NOTE: proviruses identified by geNomad were excluded here
    """
    input:
        viral_id=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_viral_contig_names.tsv"),
        contigs=os.path.join(dir["output"]["assembly"], "renamed_contigs", "{sample}_contigs_1kb.fasta")
    output:
        viral_seq=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_without_genomad_provirus.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if viral_id is not empty, extract their sequences
        if [ -s "{input.viral_id}" ]; then
            seqkit grep -i -f {input.viral_id} {input.contigs} > {output.viral_seq}
        else
            echo "viral_id file is empty"
            touch {output.viral_seq}
        fi
        """

rule extract_genomad_provirus_sequences:
    """
    get contig names of provirus identified by geNomad. TODO: double check in termina as PvVT has none
    """
    input:
        provirus_id=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_lists", "{sample}_genomad_provirus_contig_names.tsv"),
        genomad_viruses=os.path.join(dir["output"]["viral_identification"], "intermediate", "{sample}", "genomad", "{sample}_contigs_1kb_summary", "{sample}_contigs_1kb_virus.fna")
    output:
        genomad_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus.fa"),
        genomad_provirus_1kb=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus_1kb.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # if provirus_id is not empty, extract their sequences
        if [ -s "{input.provirus_id}" ]; then
            seqkit grep -i -f {input.provirus_id} {input.genomad_viruses} > {output.genomad_provirus}
            # only keep trimmed provirus with length >= 1kb 
            seqkit seq -m 1000 {output.genomad_provirus} > {output.genomad_provirus_1kb}
        else
            echo "provirus_id file is empty"
            touch {output.genomad_provirus}
            touch {output.genomad_provirus_1kb}
        fi
        """

rule combine_nonprovirus_and_genomad_provirus:
    """
    combine the non-provirus viral contigs with the genomad proviruss
    """
    input:
        genomad_provirus_1kb=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_genomad_provirus_1kb.fa"),
        viral_seq=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_without_genomad_provirus.fa")
    output:
        with_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_with_genomad_provirus.fa")
    shell:
        """
        cat {input.genomad_provirus_1kb} {input.viral_seq} > {output.with_provirus}
        """

rule checkv_identified:
    """
    run checkV on all identified viral contigs to identify additional proviruses.
    Contigs identified by VirSorter2 and Cenote-taker2 were not checked for provirus, so run checkV on the identified viral contigs and extract the trimmed provirus contigs.
    """
    input:
        with_provirus=os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "tmp", "{sample}_identified_viral_contigs_with_genomad_provirus.fa")
    output:
        checkv_dir=directory(os.path.join(dir["output"]["viral_identification"], "checkv_identified", "{sample}")),
        checkv_proviruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses.fna"), # this would NOT be a problem since CheckV produces empty provirus files
        checkv_viruses=os.path.join(dir["output"]["viral_identification"], "checkv_identified", "{sample}", "viruses.fna")
    params:
        database_path=config["checkv"]["database_path"]
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    shell:
        """
        {CHECKV}

        # if with_provirus is not empty, run checkV
        if [ -s "{input.with_provirus}" ]; then
            checkv end_to_end {input.with_provirus} {output.checkv_dir} -t {threads} -d {params.database_path}
        else
            echo "with_provirus file is empty"
            touch {output.checkv_proviruses}
            touch {output.checkv_viruses}
        fi
        """

rule combine_checkv_provirus_1kb_and_virus:
    """
    CheckV might have identified more proviruses, only keep proviruses >= 1kb.
    Combine the kept proviruses and the non-proviruses sequences
    """
    input:
        checkv_proviruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses.fna"), # this would NOT be a problem since CheckV produces empty provirus files
        checkv_viruses=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "viruses.fna")
    output:
        checkv_proviruses_1kb=os.path.join(dir["output"]["viral_identification"],"checkv_identified", "{sample}", "proviruses_1kb.fna"),
        final_viruses=os.path.join(dir["output"]["viral_identification"],"identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final.fa")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # only keep trimmed provirus with length >= 1kb 
        # deal with the possibility that no provirus identified by checkV
        if [ -s "{input.checkv_proviruses}" ]; then
            seqkit seq -m 1000 {input.checkv_proviruses} > {output.checkv_proviruses_1kb}
        else
            touch {output.checkv_proviruses_1kb}
        fi

        # combine proviruses >= 1kb with viruses
        cat {output.checkv_proviruses_1kb} {input.checkv_viruses} > {output.final_viruses}
        """

rule get_virus_sequence_length:
    """
    Get the length of final viruses
    """
    input:
        final_viruses=os.path.join(dir["output"]["viral_identification"],"identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final.fa")
    output:
        os.path.join(dir["output"]["viral_identification"], "identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final_length.txt")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    shell:
        """
        # DON'T NEED TO DEAL WITH EMPTY FASTA, seqkit WILL JUST OUTPUT HEADER
        seqkit fx2tab --length --name --header-line {input.final_viruses} > {output}
        """

rule blastn_identified:
    """blast against the bacterial genome"""
    input:
        blastdb=lambda wildcards: blastdb_index[wildcards.sample],
        contigs=os.path.join(dir["output"]["viral_identification"],"identified_viral_contig_sequences", "{sample}_identified_viral_contigs_final.fa")
    output:
        blastn=os.path.join(dir["output"]["viral_identification"], "blastn", "{sample}", "{sample}_blastn.out")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    params:
        db_path=lambda wildcards: blastdb_prefix[wildcards.sample],
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