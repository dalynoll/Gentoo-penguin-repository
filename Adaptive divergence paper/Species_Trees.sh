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


