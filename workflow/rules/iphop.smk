##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule iphop_setup:
    output:
        directory(
            config['iphop_db']
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "iphop", "iphop_setup.log"),
    resources:
        cpus=4,
    conda:
        "../envs/iphop.yaml"
    threads: 4
    shell:
        """
        iphop download --db_dir {output:q} &> {log:q}
        """


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule iphop_run:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.selected.fasta",
        ),
        iphop_db=os.path.join(
            config['iphop_db'],
            "Sept_2021_pub",
        ),
    output:
        results=directory(
            os.path.join(
                OUTPUT_FOLDER,
                "results",
                "iphop",
                "{sample}",
            ),
        )
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "iphop", "{sample}.iphop_run.log"
        ),
    resources:
        cpus=20,
    conda:
        "../envs/iphop.yaml"
    threads: 20
    shell:
        """
        iphop predict --fa_file {input.fasta:q} --out_dir {output.results} --db_dir {input.iphop_db:q} \
        --num_threads {threads} --step all &> {log:q}
        """


##########################################################################
##########################################################################