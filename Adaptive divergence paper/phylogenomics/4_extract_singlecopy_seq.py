#!/usr/bin/env python3
"""
Script to extract single-copy ortholog sequences and generate a FASTA
for each orthogroup. It is based on the file Orthogroups_SingleCopyOrthologues.txt
(produced by Orthofinder) and on the simplified proteomes (already reheadered headers)
located in a directory.
"""

import os
import sys
import csv
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Extracts single-copy sequences from orthogroups and generates a FASTA per orthogroup"
    )
    parser.add_argument(
        "-o", "--orthogroups", required=True,
        help="Path to file Orthogroups_SingleCopyOrthologues.txt (TSV format)"
    )
    parser.add_argument(
        "-f", "--fasta_dir", required=True,
        help="Directory where FASTA with simplified headers are located (eg: 1fasta_cds_curated/reheader)"
    )
    parser.add_argument(
        "-out", "--output_dir", required=True,
        help="Output directory where orthogroup FASTAs will be written (eg: 3single_copy)"
    )
    args = parser.parse_args()

    orthogroups_file = args.orthogroups
    fasta_dir = args.fasta_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # FASTA files are assumed to be named <sample>.fa
    # Example: 13A2_cds.fa, 15A2_cds.fa, etc.
    species_to_index = {}
    for file in os.listdir(fasta_dir):
        if file.endswith(".fa"):
            species = os.path.splitext(file)[0]  # extracts, for example, "13A2_cds"
            fasta_path = os.path.join(fasta_dir, file)
            print(f"Indexando {species} desde {fasta_path}...")
            species_to_index[species] = SeqIO.index(fasta_path, "fasta")

    # Opening the single-copy orthogroups (TSV) file
    with open(orthogroups_file, "r") as og_file:
        reader = csv.DictReader(og_file, delimiter="\t")
       
      
        for row in reader:                         # The header has: Orthogroup, <sample1>, <sample2>, ... according to the order used by Orthofinder.
            og_id = row["Orthogroup"]
            out_path = os.path.join(output_dir, f"{og_id}.fa")
            with open(out_path, "w") as out_f:
               
                for species in reader.fieldnames[1:]:    # For each sample (column, except the Orthogroup column)
                    gene_id_field = row.get(species, "").strip()
                    if gene_id_field == "" or gene_id_field.upper() == "NA":
                        continue
                    
                    gene_id = gene_id_field.split(",")[0].strip()  # In single-copy each cell is expected to have a unique ID.
                   
                    if species in species_to_index:     # Find the sequence in the corresponding index:
                        fasta_index = species_to_index[species]
                        if gene_id in fasta_index:
                            record = fasta_index[gene_id]
                            
                            record.id = f"{species}|{gene_id}"    # Modify the header to include the name of the sample (species)
                            record.description = ""
                            SeqIO.write(record, out_f, "fasta")
                        else:
                            print(f"Warning: {gene_id} was not found in {species}")
                    else:
                        print(f"Warning: FASTA file not found for {species}")
            print(f"Orthogroup {og_id} escrito en {out_path}")

    # Close indexes
    for idx in species_to_index.values():
        idx.close()

if __name__ == "__main__":
    main()
