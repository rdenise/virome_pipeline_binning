##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dram_setup:
    output:
        directory(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "dram_db",
            )
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "dram", "dram_setup.log"),
    resources:
        cpus=1,
    conda:
        "../envs/dram.yaml"
    threads: 1
    shell:
        """
        DRAM-setup.py prepare_databases --skip_uniref --output_dir {output:q} &> {log:q}
        DRAM-setup.py update_description_db  &>> {log:q}
        """


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dramv_annotate:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step2",
            "for-dramv",
            "final-viral-combined-for-dramv.fa",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "dram_db",
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv",
        ),
        faa=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "{sample}.faa",
        ),
    params:
        viral_affi=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "virsorter",
            "{sample}",
            "vs2-step2",
            "for-dramv",
            "viral-affi-contigs-for-dramv.tab",
        ),
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
        ),
        cutoff=cutoff_dramv,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "dramv", "{sample}.dramv_annotate.log"),
    resources:
        cpus=5,
    conda:
        "../envs/dram.yaml"
    threads: 5
    script:
        "../scripts/viral-pipeline-scripts/dram-v-wrapper.py"


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dramv_annotate_missing:
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
            "dram_db",
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "missing_annotation",
            "annotations.tsv",
        ),
        faa=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "missing_annotation",
            "{sample}.faa",
        ),
    params:
        viral_affi="",
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "missing_annotation",
        ),
        cutoff=cutoff_dramv,
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "dramv", "{sample}.dramv_annotate_missing.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/dram.yaml"
    threads: 5
    script:
        "../scripts/viral-pipeline-scripts/dram-v-wrapper.py"


##########################################################################
##########################################################################
# NOTES:
# 1. Need to think about doing the pipeline one contig by one contif or merge (as Andrey does)
# 2. In the config file or in another tabulated file have the path to all the database fasta file
# Because right now all the databases have a not similar way of being


rule dramv_distill:
    input:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "dramv",
            "annotate",
            "{sample}",
            "annotations.tsv",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "dram_db",
        ),
    output:
        amg_summary=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "dramv",
            "distill",
            "{sample}",
            "product.html",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "dramv",
            "distill",
            "{sample}",
        ),
        cutoff=cutoff_dramv,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "dramv", "{sample}.dramv_distill.log"),
    resources:
        cpus=5,
    conda:
        "../envs/dram.yaml"
    threads: 5
    shell:
        """
        [ -d {params.output_dir:q} ] && [ ! -L {params.output_dir:q} ] && rm -rf {params.output_dir:q}
        DRAM-v.py distill -i {input.tsv:q} -o {params.output_dir:q} &> {log:q}
        """


##########################################################################
##########################################################################
