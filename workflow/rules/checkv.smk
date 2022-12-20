##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_setup:
    output:
        directory(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "checkv_db",
                "checkv-db-v1.4",
            )
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "checkv", "checkv_setup.log"),
    params:
        checkv_db=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/checkv.yaml"
    threads: 1
    shell:
        """
        checkv download_database {params.checkv_db:q} &> {log:q}
        """


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_run:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step1",
            "final-viral-combined.fa",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
            "checkv-db-v1.4",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "proviruses.fna",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
            "contamination.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "virsorter",
            "{sample}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "virsorter",
            "{sample}.checkv_run.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/checkv.yaml"
    threads: 5
    shell:
        """
        checkv end_to_end {input.fasta:q} {params.output_dir:q} \
        -t {threads} -d {input.database:q} &> {log:q}

        for myfile in {output}
        do 
            touch $myfile
        done
        """


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule checkv_run_missing_annotation:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}.missing_annotation.fasta",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
            "checkv-db-v1.4",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "missing_annotation",
            "{sample}",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "missing_annotation",
            "{sample}",
            "proviruses.fna",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "missing_annotation",
            "{sample}",
            "contamination.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "checkv",
            "missing_annotation",
            "{sample}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "missing_annotation",
            "{sample}.checkv_run.log"),
    resources:
        cpus=5,
    conda:
        "../envs/checkv.yaml"
    threads: 5
    shell:
        """
        checkv end_to_end {input.fasta:q} {params.output_dir:q} \
        -t {threads} -d {input.database:q} &> {log:q}

        for myfile in {output}
        do 
            touch $myfile
        done
        """


##########################################################################
##########################################################################