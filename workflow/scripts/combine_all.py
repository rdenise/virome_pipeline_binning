from Bio import SeqIO
import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Dataframe that contains all the informations about
output_df = pd.read_table(snakemake.input.tsv_kept)

contig_reinjected = []

with open(snakemake.input.tsv_blast) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.split()[0]

        new_contig_name = f"{rstrip_line}--blast_reinjected"
        contig_reinjected.append(rstrip_line)

        output_df.at[new_contig_name, "contig_id"] = new_contig_name
        output_df.at[new_contig_name, "old_contig_id"] = rstrip_line
        output_df.at[new_contig_name, "blast_reinjected"] = "Yes"

# Fill the informations missing now the list of contigs we keep is set
dict_map_virsorter = {}

files_with_info = {
    snakemake.input.ids_virsorter_keep2_suspicious: "keep2_suspicious",
    snakemake.input.ids_virsorter_manual_check: "to_manual_check",
    snakemake.input.ids_virsorter_discarded: "discarded",
}

for file_ids in files_with_info:
    with open(file_ids) as r_file:
        r_file.readline()

        for line in r_file:
            rstrip_line = line.rstrip()
            rstrip_line, leftover = rstrip_line.split("||")

            if rstrip_line in contig_reinjected:
                new_contig_name = f"{rstrip_line}--blast_reinjected"
                dict_map_virsorter[new_contig_name] = files_with_info[file_ids]

# Fill the dataframe
list_contig2add_virsorter_cat = list(dict_map_virsorter.keys())
output_df.loc[
    output_df.contig_id.isin(list_contig2add_virsorter_cat), "virsorter_cat"
] = output_df.loc[
    output_df.contig_id.isin(list_contig2add_virsorter_cat), "contig_id"
].map(
    dict_map_virsorter
)

output_df.fillna("No", inplace=True)


with open(snakemake.output.fasta, "w") as w_file:
    SeqIO.write(SeqIO.parse(snakemake.input.fasta_kept, "fasta"), w_file, "fasta")

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.fasta_discarded

print("------------")
print(f"list only hblast reinjected: {contig_reinjected}")
print("------------")

with open(snakemake.output.fasta, "a") as w_file:
    parser = SeqIO.parse(fasta_contigs, "fasta")

    contig2new_name = output_df.set_index("old_contig_id").contig_id.to_dict()

    for contig in parser:
        contig_id = contig.id

        if contig_id in contig_reinjected :
            contig_id = contig2new_name[contig.id]

            contig.id = contig_id
            contig.name = ""
            contig.description = ""

            SeqIO.write(contig, w_file, "fasta")

output_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

###########################################################
###########################################################
