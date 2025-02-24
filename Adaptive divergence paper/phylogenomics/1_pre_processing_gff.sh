#### gff files from GCA_010090195.1 (reference genome) and GCA_030674165.1 (South Georgia gentoo assembly) were used
#### pre processing of gff to exclude pseudogenes

cat GCA_010090195.1_BGI_Ppap.V1_genomic.gff |grep -E '(gene|mRNA|exon|CDS)' > GCA_010090195.1_BGI_Ppap.V1_genomic_filtered.gff
cat GCA_030674165.1_ASM3067416v1_genomic.gff |grep -E '(gene|mRNA|exon|CDS)' > GCA_030674165.1_ASM3067416v1_genomic_filtered.gff

