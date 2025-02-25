#!/usr/bin/env python3
"""
Script to extract single-copy ortholog sequences and generate a FASTA
for each orthogroup. It uses two files:
  - Orthogroups_SingleCopyOrthologues.txt: a list (one per line) of OG IDs (no header).
  - Orthogroups.txt: full mapping file (TSV or space‐separated) where each line is:
      OGID: geneID_sample1 geneID_sample2 ... geneID_sampleN

The script requires a comma-separated list of sample names (in the same order as
the gene IDs appear in Orthogroups.txt, after the OG ID) so that it can find the corresponding
FASTA files (named <sample>.fa) in a given directory.

Excecution:
python3 3extract_single_copy_orthogroups.py \
  -sco 2orthofinder/Results_Feb24/Orthogroups/Orthogroups_SingleCopyOrthologues.txt \
  -og 2orthofinder/Results_Feb24/Orthogroups/Orthogroups.txt \
  -f 1fasta_cds_curated/reheader \
  -out 3single_copy \
  -s "13A2_cds,15A2_cds,24A1_cds,CRO05_cds,CRO07_cds,CRO08_cds,GSA10_cds,GSA7_cds,GSA9_cds,GTP06_cds,GTP07_cds,GTP12_cds,MQ9209_cds,MQ9312_cds,P0048_cds,P0052_cds,P0065_cds,P0173_cds,P0196_cds,P0200_cds,P0404_cds,P0407_cds,P0412_cds,P1885_cds,P1889_cds,P1893_cds,PAM3_cds,PAM5_cds,PAM7_cds,adel_cds,chins_cds,sgeorgia_cds"



"""

import os
import argparse

from Bio import SeqIO

def load_singlecopy_og_ids(singlecopy_path):
    """Load OG IDs from Orthogroups_SingleCopyOrthologues.txt into a set."""
    og_set = set()
    with open(singlecopy_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                og_set.add(line)
    return og_set

def parse_orthogroups(orthogroups_path, singlecopy_set):
    """
    Parse Orthogroups.txt.
    Returns a dict mapping OG_ID to a list of gene IDs (one per sample).
    Only OGs that are in the singlecopy_set are returned.
    Assumes each line is in the format:
      OGID: geneID1 geneID2 ... geneIDN
    """
    og_mapping = {}
    with open(orthogroups_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split on colon
            parts = line.split(":", 1)
            if len(parts) != 2:
                continue
            og_id = parts[0].strip()
            # Only process OGs in the singlecopy list
            if og_id not in singlecopy_set:
                continue
            # Split the rest by whitespace to get gene IDs
            gene_ids = parts[1].strip().split()
            og_mapping[og_id] = gene_ids
    return og_mapping

def main():
    parser = argparse.ArgumentParser(
        description="Extract single-copy sequences for each orthogroup and generate a FASTA per OG."
    )
    parser.add_argument(
        "-sco", "--singlecopy_file", required=True,
        help="Path to Orthogroups_SingleCopyOrthologues.txt (one OG ID per line)"
    )
    parser.add_argument(
        "-og", "--orthogroups_file", required=True,
        help="Path to Orthogroups.txt (mapping: OGID: geneID geneID ...)"
    )
    parser.add_argument(
        "-f", "--fasta_dir", required=True,
        help="Directory where FASTA files with simplified headers are located (e.g., 1fasta_cds_curated/reheader)"
    )
    parser.add_argument(
        "-out", "--output_dir", required=True,
        help="Output directory where OG FASTAs will be written (e.g., 3single_copy)"
    )
    parser.add_argument(
        "-s", "--samples", required=True,
        help=("Comma-separated list of sample names, in the order corresponding to the columns in Orthogroups.txt "
              "(excluding the OG ID). For example: "
              "\"13A2_cds,15A2_cds,24A1_cds,CRO05_cds,...,sgeorgia_cds\"")
    )
    args = parser.parse_args()

    singlecopy_file = args.singlecopy_file
    orthogroups_file = args.orthogroups_file
    fasta_dir = args.fasta_dir
    output_dir = args.output_dir
    sample_list = [s.strip() for s in args.samples.split(",")]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load single-copy OG IDs into a set
    singlecopy_set = load_singlecopy_og_ids(singlecopy_file)
    print(f"Loaded {len(singlecopy_set)} single-copy OG IDs.")

    # Parse the full Orthogroups.txt mapping (only keeping OGs in singlecopy_set)
    og_mapping = parse_orthogroups(orthogroups_file, singlecopy_set)
    print(f"Found {len(og_mapping)} OG mappings in Orthogroups.txt for single-copy OGs.")

    # Index FASTA files for each sample (assumed to be named <sample>.fa)
    species_to_index = {}
    for file in os.listdir(fasta_dir):
        if file.endswith(".fa"):
            species = os.path.splitext(file)[0]  # e.g., "13A2_cds"
            fasta_path = os.path.join(fasta_dir, file)
            print(f"Indexing {species} from {fasta_path}...")
            species_to_index[species] = SeqIO.index(fasta_path, "fasta")

    # Process each orthogroup
    for og_id, gene_ids in og_mapping.items():
        # Check that the number of gene IDs matches the number of samples
        if len(gene_ids) != len(sample_list):
            print(f"Warning: In orthogroup {og_id}, number of gene IDs ({len(gene_ids)}) does not match sample list length ({len(sample_list)}). Skipping this OG.")
            continue
        out_path = os.path.join(output_dir, f"{og_id}.fa")
        with open(out_path, "w") as out_f:
            # For each sample, get the corresponding gene ID (by column order)
            for i, sample in enumerate(sample_list):
                gene_id = gene_ids[i].strip()
                # If gene_id is empty or "NA", skip it
                if gene_id == "" or gene_id.upper() == "NA":
                    continue
                if sample in species_to_index:
                    fasta_index = species_to_index[sample]
                    if gene_id in fasta_index:
                        record = fasta_index[gene_id]
                        # Modify header to include sample name and gene id
                        record.id = sample.split("_cds")[0]
                        record.description = ""
                        SeqIO.write(record, out_f, "fasta")
                    else:
                        print(f"Warning: {gene_id} not found in {sample}")
                else:
                    print(f"Warning: FASTA file not found for sample {sample}")
        print(f"Orthogroup {og_id} written to {out_path}")

    # Close all FASTA indexes
    for idx in species_to_index.values():
        idx.close()

if __name__ == "__main__":
    main()
