##########################################################
######################## Run AGAT ########################
##########################################################

### Extract CDS from gff file:
file="sample1
sample2
sample3"
for i in $file
do
agat_sp_extract_sequences.pl -g $anotation/genomic.gff -f $genomes/${i}_maskedPygo.fa -o $cds/1_agat/${i}_cds.fa -t cds
done
##### gff from GCA_010090195.1 (reference genome) and GCA_030674165.1 (South Georgia gentoo assembly) were used


#####################################################################################
######################## Modify header of agat extracted CDS ########################
#####################################################################################

### rename header of gene name creating a id2gene.map
awk -F '\t' '
  $3 == "CDS" {
    match($9, /Parent=([^;]+)/, par)
    match($9, /gene=([^;]+)/,   gen)
    if (par[1] != "" && gen[1] != "") {
      # Imprimimos Parent [tab] gene
      print par[1] "\t" gen[1]
    }
  }
' genomic.gff > id2gene.map


### Create rename_fasta.awk:
BEGIN {
  FS = "\t"
}
# Read map ID -> gene in the phase fase NR==FNR
NR == FNR {
  dic[$1] = $2
  next
}
# Then, when process the FASTA file:
{
  if ($0 ~ /^>/) {
    # Remove '>' 
    header = substr($0, 2)
    # The real ID is usually up to the first space
    split(header, arr, " ")
    oldID = arr[1]

    if (oldID in dic) {
      # If it is in the dictionary, we use the gene name
      print ">" dic[oldID]
    } else {
      # If not, we leave the header as is
      print ">" oldID
    }
  } else {
        print
  }
}

### Run:
awk -f rename_fasta.awk id2gene.map 1_agat/${i}_cds.fa > 2_cds_reheader/${i}_cds_renamed.fa



#################################################################################################################
##################### Search common gen names with South Georgea assembly (GCA_030674165.1) #####################
#################################################################################################################
cat *.fa | grep '^>' | tr '[:lower:]' '[:upper:]' | sort | uniq -c | grep ' 67 ' | cut -d'>' -f2 > core_genes_corrected.txt
### 10,299 genes
### 67 samples: 65 gentoo penguins (including 2 from Macquarie and 1 from South Georgias), one chinstrap and one adelie penguin

##### Extract common genes:
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



############################################################################
######################## Align sequences with MAFFT ########################
############################################################################

#!/bin/bash
mkdir -p 4_align_mafft

for fasta in 3_seq_genes/*.fa; do
    gene=$(basename "$fasta" .fa)
    echo "Alineando $gene con MAFFT..."

    mafft --auto "$fasta" > "4_align_mafft/${gene}_aligned.fa"
done





