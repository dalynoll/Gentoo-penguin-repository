#!/usr/bin/env python3

import os
indir = "2_cds_reheader"
outdir = "3_seq_genes"
os.makedirs(outdir, exist_ok=True)

core_file = "3_seq_genes/core_genes_corrected.txt"
with open(core_file, "r") as cf:
    core_set = { line.strip().upper() for line in cf if line.strip() }

for filename in os.listdir(indir):
    if filename.endswith("_cds_renamed.fa"):
        sample = filename.replace("_cds_renamed.fa", "")
        print(f"Processing {filename}...")

        with open(os.path.join(indir, filename), "r") as infile:
            current_gene = None
            seq_lines = []
            
            for line in infile:
                line = line.rstrip("\n")
                if line.startswith(">"):
 
                    if current_gene and current_gene.upper() in core_set and seq_lines:
                        out_path = os.path.join(outdir, f"{current_gene.upper()}.fa")
                        with open(out_path, "a") as gene_file:
                            gene_file.write(f">{sample}\n")
                            gene_file.write("\n".join(seq_lines) + "\n")

                    current_gene = line[1:].strip()  # remover '>'
                    seq_lines = []
                else:
                    seq_lines.append(line)
    
            if current_gene and current_gene.upper() in core_set and seq_lines:
                out_path = os.path.join(outdir, f"{current_gene.upper()}.fa")
                with open(out_path, "a") as gene_file:
                    gene_file.write(f">{sample}\n")
                    gene_file.write("\n".join(seq_lines) + "\n")

