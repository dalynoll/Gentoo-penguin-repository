
# **Workflow VCF analysis gentoo penguin population structure**

## 1 Outlier SNPs exclusion
vcftools --gzvcf $vcf_64 --bed final_positions_noOutliers.bed --recode --out vcf_64_noOutliers
###### IMPORTANT: The BED file contains SNPs located within RAiSD outlier windows

## 2 Minor allele frequency and Pruning by LD

plink --vcf $vcf_64_noOutliers --maf 0.01 --indep-pairwise 50 10 0.1 --out gentoo64maf001 --set-missing-var-ids @:# --double-id --allow-extra-chr

plink --vcf $vcf_64_noOutliers --extract gentoo64maf001.prune.in --double-id  --allow-extra-chr --set-missing-var-ids @:# --recode vcf --out gentoo64_noOutliers_maf001_50_10_01



## 3 PCA

plink --file $gentoo64_noOutliers_maf001_50_10_01 --out pca_gentoo_64_noOutliers --pca --allow-extra-chr


## 4 Run ADMIXTURE


#### ADMIXTURE input file 
plink --vcf $gentoo64_noOutliers_maf001_50_10_01 --extract gentoo64maf001.prune.in --double-id --allow-extra-chr --set-missing-var-ids @:# --recode 12 --out gentoo64_MAF001_50_10_01


```bash
for K in 1 2 3 4 5 6 7 8 9 10 
do 
	admixture --cv gentoo64_noOutliers_maf001_50_10_01.ped $K -j8 | tee log${K}.out 
done 
```

## 5. Heterozygosity
vcftools --vcf $gentoo64_noOutliers_maf001_50_10_01 --het --out gentoo64_noOut_het


## 6. FST - DXY using pixy

###### First, a new variant calling of variant and invariant sites was performed:

cat chromosome.list.selection | xargs -I {} -n 1 -P 40 sh -c 'bcftools mpileup -Ou -f GCA_010090195.1_BGI_Ppap.V1_genomic.fna -b bam_list -d 5000 -q 10 -Q 20 -a SP,DP --skip-indels -r {} --rf 2 | bcftools call -f GQ -m -Oz -o {}.vcf.gz'


#### windows 5kb (genome divergence Figure 1f)
pixy --stats pi dxy fst --vcf $VCF_var_invar --populations pop_list.txt --window_size 5000 --n_cores 20 --fst_type wc 


#### windows 200kb (detection of genomic areas with high relative (FST) and low absolute (DXY) differentiation, Figure 3)

pixy --stats pi dxy fst --vcf $VCF_var_invar --populations pop_list.txt --window_size 200000 --n_cores 20 --fst_type wc 

