##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being
# Important:: Need to maybe separate the contig inside the BIG contig file because each could be
# a bacteria or a virus by itself


rule vcontact2_preprocess:
    input:
        proteins_fasta=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "dramv",
                "annotate",
                "{sample}",
                "merge",
                "{sample}.faa",
            ),
            sample=CONTIGS_DICT.keys(),
        ),
        database=os.path.join(
            DB_DICT["protein"]["ICTV"]["path"],
            DB_DICT["protein"]["ICTV"]["file"],
        ) if "ICTV" in DB_DICT["protein"] else False,
        db_genomes=os.path.join(
            DB_DICT["fasta"]["ICTV"]["path"],
            DB_DICT["fasta"]["ICTV"]["file"],
        ) if "ICTV" in DB_DICT["fasta"] else False,
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "query",
            "merge.rename.faa",
        ),
        csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "query",
            "merge.proteins.csv",
        ),
        fasta_low=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "query",
            "merge.rename.low.faa",
        ),
    params:
        fasta_contig=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "viral_contigs",
                "{sample}.selected.fasta",
            ), 
        sample=CONTIGS_DICT.keys(),
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vcontact2",
            "contig.vcontact2_preprocess.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/protein_csv.py"


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule vcontact2:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "query",
            "merge.rename.faa",
        ),
        protein_csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "query",
            "merge.proteins.csv",
        ),
    output:
        csv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
            "genome_by_genome_overview.csv",
        ),
    params:
        vcontact2_db="ProkaryoticViralRefSeq211-Merged",
        rel_mode="Diamond",
        pcs_mode="MCL",
        vcs_mode="ClusterONE",
        c1_bin="clusterone",
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "vcontact2",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "vcontact2", "vcontact2.log"),
    resources:
        cpus=25,
    conda:
        "../envs/vcontact2.yaml"
    threads: 25
    shell:
        """
        vcontact2 --raw-proteins {input.fasta:q} --rel-mode {params.rel_mode} \
        --proteins-fp {input.protein_csv:q} --db {params.vcontact2_db} \
        --pcs-mode {params.pcs_mode} --vcs-mode {params.vcs_mode} \
        --c1-bin {params.c1_bin:q} \
        --output-dir {params.output_dir:q} -t {threads} &> {log:q}
        """


##########################################################################
##########################################################################
