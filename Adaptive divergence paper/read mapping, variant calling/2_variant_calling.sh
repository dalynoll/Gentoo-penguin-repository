######################################
########### Variant calling ##########
######################################

cat chromosome.list.selection | xargs -I {} -n 1 -P 40 sh -c 'bcftools mpileup -Ou \
-f GCA_010090195.1_BGI_Ppap.V1_genomic.fna -b bam_list -d 5000 -q 10 -Q 20 -a SP,DP \
--skip-indels -r {} --rf 2 | bcftools call -f GQ -vm -Oz -o {}.vcf.gz'



### concatenate scaffolds vcf
bcftools concat -f concat.list --threads 10 -Ov -o unfiltered_gentoo.vcf
