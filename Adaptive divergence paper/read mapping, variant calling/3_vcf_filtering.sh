  # # # # # ## # # # # ## # # # # # # # # #
# # # # # # Variant Normalization # # # # # #
  # # # # # ## # # # # ## # # # # # # # # #

bcftools norm --check-ref w -f GCA_010090195.1_BGI_Ppap.V1_genomic.fna \
-o unfiltered_gentoo_norm.vcf.gz -Oz --threads 8  unfiltered_gentoo.vcf


  # # # # # # # # # # ## # # # # # # # # # # # # #
# # # # # # Filtering of SNPs and Indels # # # # # #
  # # # # # # # # # # ## # # # # # # # # # # # # #
bcftools filter --SnpGap 5 --IndelGap 5 -Oz -o gentoo_indels.vcf.gz  unfiltered_gentoo_norm.vcf

