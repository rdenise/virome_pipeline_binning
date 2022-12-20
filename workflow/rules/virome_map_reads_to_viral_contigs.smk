##########################################################################
##########################################################################

rule build_index:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "viral_contigs",
            "{sample}_contigs.selected.fasta",
        ),
    output:
        index=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "bowtie2",
            "{sample}",
            "index",
            "{sample}_contigs.selected.fasta.1.bt2",
        )
    params:
        index=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "bowtie2",
            "{sample}",
            "index",
            "{sample}_contigs.selected.fasta",
        ),     
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2_build",
            "{sample}.bowtie2_build.log",
        ),
    resources:
        cpus=1,
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {input.fasta:q} {params.index:q} &> {log:q}
        """

##########################################################################
##########################################################################


rule map_reads:
    input:
        index=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "bowtie2",
            "{sample}",
            "index",
            "{sample}_contigs.selected.fasta.1.bt2",
        ),
        r1=os.path.join(
            FASTQ_FOLDER,
            "{sample}_1.fastq.gz"
        ),
        r2=os.path.join(
            FASTQ_FOLDER,
            "{sample}_2.fastq.gz"
        ),
    output:
        temp(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "bowtie2",
                "{sample}",                
                "{sample}.sam",
            )
        )
    params:
        index=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "bowtie2",
            "{sample}",
            "index",
            "{sample}_contigs.selected.fasta",
        ),        
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2",
            "{sample}.bowtie2.log",
        ),
    conda:
        "../envs/bowtie2.yaml"
    threads:
        10
    shell:
        """
        bowtie2 --end-to-end -x {params.index:q} -1 {input.r1:q} -2 {input.r2:q} -p {threads} -S {output:q} &> {log:q}
        """


##########################################################################
##########################################################################

rule samtools_view:
    input:
        sam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "bowtie2",
            "{sample}",
            "{sample}.sam",
        )
    output:
        bam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}.sam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.view.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -S -b -@ {threads} -o {output.bam:q} {input.sam:q} &> {log:q}
        """

##########################################################################
##########################################################################
rule samtools_sort:
    input:
        bam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}.sam",
        ),
    output:
        sortedbam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.sort.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -m 20G -o {output.sortedbam:q} {input.bam:q} &> {log:q}
        """

##########################################################################
##########################################################################


rule samtools_index:
    input:
        sortedbam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        idxstats=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.idxstats",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.idxstats.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index -@ {threads} -b {input.sortedbam:q} &> {log:q}
        samtools idxstats {input.sortedbam:q} > {output.idxstats:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule samtools_flagstat:
    input:
        sortedbam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        flags=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.flags",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.flagstat.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat -@ {threads} {input.sortedbam:q} > {output.flags:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule samtools_pileup:
    input:
        sortedbam=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        pileup=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.pileup",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.mpileup.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools mpileup {input.sortedbam:q} > {output.pileup:q} 2> {log:q}
        """

##########################################################################
##########################################################################