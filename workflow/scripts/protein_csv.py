from Bio import SeqIO
import sys, os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

with open(snakemake.output.csv, "wt") as w_file:
    with open(snakemake.output.fasta, "wt") as fasta_file:
        with open(snakemake.output.fasta_low, "wt") as fasta_file_low:

            header = "contig_id,protein_id,keywords"
            w_file.write(f"{header}\n")

            # ICTV
            if snakemake.input.database:
                parser = SeqIO.parse(snakemake.input.database, "fasta")

                contig_db_parser = SeqIO.parse(snakemake.input.db_genomes, "fasta")
                dict_contig_db_len = {contig.id: len(contig.seq) for contig in contig_db_parser}

                # Parse database and remove case only if I can get the contig length
                for protein in parser:
                    protein_id = protein.id
                    contig_id = protein.description.split(" ")[8]
                    keyword = "ICTV"

                    # print(protein.description.split(" "))

                    if dict_contig_db_len[contig_id] >= 10000:
                        w_file.write(f"{contig_id},{protein_id},{keyword}\n")

                        protein.name = ""
                        protein.description = ""

                        SeqIO.write(protein, fasta_file, "fasta")
                    else:
                        protein.name = ""
                        protein.description = ""

                        SeqIO.write(protein, fasta_file_low, "fasta")
                                           
            num_files = len(snakemake.input.proteins_fasta)
            for index in range(num_files):
                protein_file = snakemake.input.proteins_fasta[index]
                parser = SeqIO.parse(protein_file, "fasta")

                contig_parser = SeqIO.parse(snakemake.params.fasta_contig[index], "fasta")

                dict_contig_len = {contig.id: len(contig.seq) for contig in contig_parser}

                print(snakemake.params.fasta_contig[index])
                for protein in parser:
                    protein_id = protein.id
                    contig_id = "_".join(protein_id.split("_")[:-1])
                    keyword = "None"

                    if dict_contig_len[contig_id] >= 10000:
                        w_file.write(f"{contig_id},{protein_id},{keyword}\n")

                        protein.name = ""
                        protein.description = ""

                        SeqIO.write(protein, fasta_file, "fasta")
                    else:
                        protein.name = ""
                        protein.description = ""

                        SeqIO.write(protein, fasta_file_low, "fasta")

###########################################################
###########################################################
