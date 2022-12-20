# Module containing all the ncbi-blast related rules


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# Advantage of the ocncatenation, no nee to do it after by yourself
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule blastn:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER, "databases", "viral_contigs", "{sample}.selected.fasta"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "virus",
            "{sample}",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.txt",
        ),
    params:
        database=lambda wildcards: os.path.join(
            DB_DICT["fasta"][wildcards.database]["path"],
            DB_DICT["fasta"][wildcards.database]["file"],
        ),
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "human", "{sample}_tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        options_blast="-num_alignments 25000",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.log",
        ),   
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 5
    script:
        "../scripts/blastn_wrapper.py"


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# Advantage of the ocncatenation, no nee to do it after by yourself
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule blastn_preprocess:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER, "processing_files", "discarded_contig", "{sample}.discarded.fasta"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "discarded",
            "{sample}",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.txt",
        ),
    params:
        database=lambda wildcards: os.path.join(
            DB_DICT["discarded"][wildcards.database]["path"],
            DB_DICT["discarded"][wildcards.database]["file"],
        ),
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "human", "{sample}_tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        options_blast="-num_alignments 25000",        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "discarded",
            "{sample}.evalue_{evalue}.{database}.blastn.outfmt6.log",
        ),
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 5
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
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "{sample}.nt.human.blastn.outfmt6.log",
        ),
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 20
    script:
        "../scripts/blastn_wrapper.py"
