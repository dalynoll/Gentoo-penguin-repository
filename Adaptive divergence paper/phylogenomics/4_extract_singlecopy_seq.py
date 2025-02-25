#!/usr/bin/env python3
"""
Script to extract single-copy ortholog sequences and generate a FASTA
for each orthogroup. It is based on the file Orthogroups_SingleCopyOrthologues.txt
(produced by Orthofinder) and on the simplified proteomes (already reheadered headers)
located in a directory.

Excecution:
python3 3extract_single_copy_orthogroups.py \
  -o 2orthofinder/Results_Feb24/Orthogroups/Orthogroups_SingleCopyOrthologues.txt \
  -f 1fasta_cds_curated/reheader \
  -out 3single_copy \
  -s "13A2_cds,15A2_cds,24A1_cds,CRO05_cds,CRO07_cds,CRO08_cds,GSA10_cds,GSA7_cds,GSA9_cds,GTP06_cds,GTP07_cds,GTP12_cds,MQ9209_cds,MQ9312_cds,P0048_cds,P0052_cds,P0065_cds,P0173_cds,P0196_cds,P0200_cds,P0404_cds,P0407_cds,P0412_cds,P1885_cds,P1889_cds,P1893_cds,PAM3_cds,PAM5_cds,PAM7_cds,adel_cds,chins_cds,sgeorgia_cds"

"""

import os
import csv
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Extracts single-copy sequences from orthogroups and generates a FASTA per orthogroup"
    )
    parser.add_argument(
        "-o", "--orthogroups", required=True,
        help="Path to file Orthogroups_SingleCopyOrthologues.txt (TSV format, no header)"
    )
    parser.add_argument(
        "-f", "--fasta_dir", required=True,
        help="Directory where FASTA with simplified headers are located (e.g., 1fasta_cds_curated/reheader)"
    )
    parser.add_argument(
        "-out", "--output_dir", required=True,
        help="Output directory where orthogroup FASTAs will be written (e.g., 3single_copy)"
    )
    parser.add_argument(
        "-s", "--samples", required=True,
        help=("Comma-separated list of sample names, in the order corresponding to the columns in the orthogroups file "
              "(excluding the first column which is the orthogroup ID).")
    )
    args = parser.parse_args()

    orthogroups_file = args.orthogroups
    fasta_dir = args.fasta_dir
    output_dir = args.output_dir
    sample_list = [s.strip() for s in args.samples.split(",")]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # FASTA files are assumed to be named <sample>.fa
    # Example: 13A2_cds.fa, 15A2_cds.fa, etc.
    species_to_index = {}
    for file in os.listdir(fasta_dir):
        if file.endswith(".fa"):
            species = os.path.splitext(file)[0]  # extracts, for example, "13A2_cds"
            fasta_path = os.path.join(fasta_dir, file)
            print(f"Indexing {species} from {fasta_path}...")
            species_to_index[species] = SeqIO.index(fasta_path, "fasta")

    # Process the orthogroups file (no header, TSV)
    with open(orthogroups_file, "r") as og_file:
        reader = csv.reader(og_file, delimiter="\t")
        for row in reader:
            if not row:
                continue
            # The first column is the orthogroup ID
            og_id = row[0]
            out_path = os.path.join(output_dir, f"{og_id}.fa")
            with open(out_path, "w") as out_f:
                # Ensure that the number of columns (excluding the first) matches the sample list length
                if len(row) - 1 != len(sample_list):
                    print(f"Warning: In orthogroup {og_id}, number of columns ({len(row)-1}) does not match sample list length ({len(sample_list)}). Skipping this OG.")
                    continue
                # For each sample, in order
                for i, sample in enumerate(sample_list):
                    gene_id_field = row[i+1].strip()
                    if gene_id_field == "" or gene_id_field.upper() == "NA":
                        continue
                    # In single-copy, each cell is expected to have a unique ID; if there are more, take the first.
                    gene_id = gene_id_field.split(",")[0].strip()
                    if sample in species_to_index:
                        fasta_index = species_to_index[sample]
                        if gene_id in fasta_index:
                            record = fasta_index[gene_id]
                            # Modify header: include sample name and gene ID
                            record.id = f"{sample}|{gene_id}"
                            record.description = ""
                            SeqIO.write(record, out_f, "fasta")
                        else:
                            print(f"Warning: {gene_id} not found in {sample}")
                    else:
                        print(f"Warning: FASTA file not found for sample {sample}")
            print(f"Orthogroup {og_id} written to {out_path}")

    # Close all indexes
    for idx in species_to_index.values():
        idx.close()

if __name__ == "__main__":
    main()
