import sys
from Bio import SeqIO
import pandas as pd
import re

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

tsv_missing = snakemake.input.tsv_missing

faa_missing = snakemake.input.faa_missing

tsv_virsorter = snakemake.input.tsv_virsorter

faa_virsorter = snakemake.input.faa_virsorter

transl_table = snakemake.input.tsv
transl_table = pd.read_table(transl_table)

##########################################################################

missing_df = pd.read_table(tsv_missing, index_col=0)
virsorter_df = pd.read_table(tsv_virsorter, index_col=0)

virsorter_df.index.name = 'prot_id'
virsorter_df = virsorter_df.reset_index().replace({"__":"||", "-cat_[0-9]+":""}, regex=True)

dict_translate = transl_table.set_index("id_virsorter").contig_id.to_dict()
dict_translate = {re.escape(k):v for k, v in dict_translate.items()}
 
virsorter_df = virsorter_df.replace(dict_translate, regex=True)

virsorter_df = virsorter_df[virsorter_df.scaffold.isin([*dict_translate.values()])].set_index("prot_id")
virsorter_df.index.name = ''

concat_df = pd.concat([virsorter_df, missing_df])
concat_df.to_csv(snakemake.output.tsv, sep="\t")

all_contigs = concat_df.scaffold.tolist()

##########################################################################

with open(snakemake.output.fasta, "wt") as w_file:
    for faa_file in [faa_virsorter, faa_missing]:
        parser = SeqIO.parse(faa_file, "fasta")
        for protein in parser:
            contig_name = protein.id.split("-cat")[0].replace("__", "--")
            if contig_name in all_contigs:
                # Need to change only the virsorter genes
                if "partial" in protein.id or "full" in protein.id or "lt2gene" in protein.id:
                    protein_id = re.sub("-cat_[0-9]+","",protein.id.replace("__", "--"))
                    
                    print(f"old name: {protein.id}")
                    protein.id = protein_id
                    print(f"new name {protein.id}")
                    print("---------------------")

                protein.name = protein.description = ""
                SeqIO.write(protein, w_file, "fasta")
