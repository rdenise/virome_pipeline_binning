from Bio import SeqIO
import sys
import pandas as pd

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Get all the names from the deepvirfinder list and remove redondant name
ids_kept = snakemake.input.tsv

ids_kept = pd.read_table(ids_kept).old_contig_id.tolist()

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contig


with open(snakemake.output.fasta, "w") as w_file:
    parser_fasta = SeqIO.parse(fasta_contigs, "fasta")

    for seq in parser_fasta:
        if seq.id not in ids_kept:
            SeqIO.write(seq, w_file, "fasta")

###########################################################
###########################################################
