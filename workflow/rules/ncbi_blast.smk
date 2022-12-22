# Module containing all the ncbi-blast related rules


##########################################################################
##########################################################################


rule blastn:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "{sample}.filtered.sorted.fasta",
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "preprocess",
            "{sample}",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.txt",
        ),
    params:
        database=lambda wildcards: os.path.join(
            DB_DICT["fasta"][wildcards.database]["path"],
            DB_DICT["fasta"][wildcards.database]["file"],
        ),
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "preprocess", "{sample}_tmp_{database}"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        options_blast="-num_alignments 25000",
        max_len=2000000,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "preprocess",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.log",
        ),   
    resources:
        cpus=20,
    conda:
        "../envs/blast.yaml"
    threads: 20
    script:
        "../scripts/blastn_wrapper.py"


##########################################################################
##########################################################################


rule blastn_human:
    input:
        contig=lambda wildcards: os.path.join(
            CONTIGS_FOLDER,
            CONTIGS_DICT[wildcards.sample]["file"],
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "human",
            "{sample}.nt.human.blastn.outfmt6.txt",
        ),
    params:
        database=blast_database,
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "human", "{sample}_tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        options_blast=(
                    "-word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -taxids 9606 "
                    "-min_raw_gapped_score 100 -perc_identity 90 -soft_masking true -max_target_seqs 10 "
        ),
        max_len=False,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "{sample}.nt.human.blastn.outfmt6.log",
        ),
    resources:
        cpus=20,
    conda:
        "../envs/blast.yaml"
    threads: 20
    script:
        "../scripts/blastn_wrapper.py"
