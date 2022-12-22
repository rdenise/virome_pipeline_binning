import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Get the informations from virsorter
deepvirfinder_file = snakemake.input.txt
deepvirfinder_file_df = pd.read_table(deepvirfinder_file)

# Good hit virfinder (form Camarillo-Gurerro_2021_Cell): score > 0.9 and pvalue < 0.01
deepvirfinder_file_df[
    (deepvirfinder_file_df.score > 0.9) & (deepvirfinder_file_df.pvalue < 0.01)
].name.to_csv(snakemake.output.txt, index=False)

###########################################################
###########################################################
