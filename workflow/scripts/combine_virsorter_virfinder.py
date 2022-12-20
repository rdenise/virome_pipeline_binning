from Bio import SeqIO
import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Get the translation step-1 (column 9) step-2 (seqname) of virsorter2
step2_step1_translate = pd.read_table(snakemake.input.final_score, index_col=9).seqname.to_dict()

# List that will contains all the contigs to filter
all_contig_ids = []

# Dataframe that contains all the informations about
output_df = pd.DataFrame(columns=["contig_id", "old_contig_id", "id_virsorter", "virsorter_cat", "deepvirfinder"])

# Get all the names from the virsorter keep2 list
ids_virsorter_keep2 = snakemake.input.ids_virsorter_keep2_checked

with open(ids_virsorter_keep2) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()
        # Because in leftover there will be the number 0 or 1 or more depending on if more than one prophage on contig 
        rstrip_line, leftover = rstrip_line.split("||")

        # Here the name from the file is the name of step 1 not step 2 of virsorter so need to translate before changing || to --
        new_contig_name = step2_step1_translate[line.rstrip()].replace("||", "--")

        if rstrip_line not in all_contig_ids:
            all_contig_ids.append(rstrip_line)

        output_df.at[new_contig_name, "contig_id"] = new_contig_name
        output_df.at[new_contig_name, "old_contig_id"] = rstrip_line
        # Here we want the name of step 2 not step 1 of virsorter
        output_df.at[new_contig_name, "id_virsorter"] = step2_step1_translate[line.rstrip()]
        output_df.at[new_contig_name, "virsorter_cat"] = "keep2_checked"

# Get all the names from the virsorter keep1 list and remove redondant name
ids_virsorter_keep1 = snakemake.input.ids_virsorter_keep1

with open(ids_virsorter_keep1) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()
        rstrip_line, leftover = rstrip_line.split("||")

        # Here the name from the file is the name of step 1 not step 2 of virsorter so need to translate before changing || to --
        new_contig_name = step2_step1_translate[line.rstrip()].replace("||", "--")

        if rstrip_line not in all_contig_ids:
            all_contig_ids.append(rstrip_line)

            output_df.at[new_contig_name, "contig_id"] = new_contig_name
            output_df.at[new_contig_name, "old_contig_id"] = rstrip_line
            output_df.at[new_contig_name, "id_virsorter"] = step2_step1_translate[line.rstrip()]
            output_df.at[new_contig_name, "virsorter_cat"] = "keep1"

# Get all the names from the deepvirfinder list and remove redondant name
ids_virfinder = snakemake.input.ids_virfinder
contig_virfinder_only = []

with open(ids_virfinder) as r_file:
    r_file.readline()

    for line in r_file:
        rstrip_line = line.rstrip()

        if rstrip_line not in all_contig_ids:
            new_contig_name = f"{rstrip_line}--deepvirfinder-only"
            contig_virfinder_only.append(rstrip_line)

            output_df.at[new_contig_name, "contig_id"] = new_contig_name
            output_df.at[new_contig_name, "old_contig_id"] = rstrip_line
            output_df.at[new_contig_name, "deepvirfinder"] = "Yes"
        else: 
            output_df.loc[output_df.old_contig_id == rstrip_line, "deepvirfinder"] = "Yes"


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

            if rstrip_line in contig_virfinder_only:
                new_contig_name = f"{rstrip_line}--deepvirfinder-only"
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

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contigs_virsorter

with open(snakemake.output.fasta, "w") as w_file:
    index_db = SeqIO.index(fasta_contigs, "fasta")

    idvir2newname = output_df[output_df.virsorter_cat.isin(["keep2_checked", "keep1"])].set_index('id_virsorter').contig_id.to_dict()

    for contig_name in idvir2newname:
        contig = index_db[contig_name]

        contig.id = idvir2newname[contig.id]
        contig.name = ""
        contig.description = ""

        SeqIO.write(contig, w_file, "fasta")

# Parse the fasta of the contig and create the new one
fasta_contigs = snakemake.input.contigs_deepvirfinder

print("------------")
print(f"list only deepvirfinder: {contig_virfinder_only}")
print("------------")

with open(snakemake.output.fasta, "a") as w_file:
    parser = SeqIO.parse(fasta_contigs, "fasta")

    contig2new_name = output_df.set_index("old_contig_id").contig_id.to_dict()

    for contig in parser:
        contig_id = contig.id

        if contig_id in contig_virfinder_only :
            contig_id = contig2new_name[contig.id]

            contig.id = contig_id
            contig.name = ""
            contig.description = ""

            SeqIO.write(contig, w_file, "fasta")

output_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

###########################################################
###########################################################
