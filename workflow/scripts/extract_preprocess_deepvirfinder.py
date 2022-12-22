from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Dataframe that contains all the informations about
output_df = pd.read_table(snakemake.input.tsv)

contig_reinjected = []

with open(snakemake.input.tsv_blast) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.split()[0]

        new_contig_name = f"{rstrip_line}--deepvirfinder_detected"
        contig_reinjected.append(rstrip_line)

        output_df.at[new_contig_name, "contig_id"] = new_contig_name
        output_df.at[new_contig_name, "old_contig_id"] = rstrip_line

with open(snakemake.output.fasta, "w") as w_file:
    SeqIO.write(SeqIO.parse(snakemake.input.fasta_kept, "fasta"), w_file, "fasta")

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.fasta

with open(snakemake.output.fasta, "w") as w_file_good:
    with open(snakemake.output.discarded, "w") as w_file_bad:
        parser = SeqIO.parse(fasta_contigs, "fasta")

        contig2new_name = output_df.set_index("old_contig_id").contig_id.to_dict()

        for contig in parser:
            contig_id = contig.id

            if contig_id in contig_reinjected :
                contig_id = contig2new_name[contig.id]

                contig.id = contig_id
                contig.name = ""
                contig.description = ""

                SeqIO.write(contig, w_file_good, "fasta")
            
            else:
                SeqIO.write(contig, w_file_bad, "fasta")

###########################################################
###########################################################
