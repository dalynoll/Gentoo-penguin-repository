
# **Workflow VCF analysis gentoo penguin resequencing**

## 1 Outlier SNPs exclusion
vcftools --gzvcf $vcf_64 --bed final_positions_noOutliers.bed --recode --out vcf_64_noOutliers

## 2 Minor allele frequency and Pruning by LD

plink --vcf $vcf_64_noOutliers --maf 0.01 --indep-pairwise 50 10 0.1 --out gentoo64maf001 --set-missing-var-ids @:# --double-id --allow-extra-chr

plink --vcf $vcf_64_noOutliers --extract gentoo64maf001.prune.in --double-id  --allow-extra-chr --set-missing-var-ids @:# --recode vcf --out gentoo64maf001_50_10_01


#### to use in ADMIXTURE
plink --vcf $vcf_64_noOutliers --extract gentoo64maf001.prune.in --double-id --allow-extra-chr --set-missing-var-ids @:# --recode 12 --out gentoo64_MAF001_50_10_01


## 3 PCA

plink --file $vcf_64_noOutliers --out pca_gentoo_64_noOutliers --pca --allow-extra-chr


## 4 Run ADMIXTURE


for K in 1 2 3 4 5 6 7 8 9 10
do
time admixture --cv vcf_64_noOutliers.ped $K -j8 | tee log${K}.out
done


## 5. Heterozygosity
vcftools --vcf $vcf_64_noOutliers --het --out gentoo64_noOut_het


## 6. FST - DXY using pixy
#### windows 5kb (genome divergence Figure 1f)
pixy --stats pi dxy fst --vcf $VCF_var_invar --populations pop_list.txt --window_size 5000 --n_cores 20 --fst_type wc 


#### windows 200kb (detection of genomic areas with high relative (FST) and low absolute (DXY) differentiation, Figure 3)

pixy --stats pi dxy fst --vcf $VCF_var_invar --populations pop_list.txt --window_size 200000 --n_cores 20 --fst_type wc 


## 7. RAiSD
### Requirements: separate the vcf by linbeage and scaffolds/chromosome. The runs are independent for each scaffold/chromosome

#### for each lineage, do:

###### split VCF by scaffold/chromosome
list=$(cat "scaff_name.txt")

for i in $list;
do vcftools  --vcf  $VCF_FILE  --chr $i  --recode --recode-INFO-all --out  $i&
done
wait


###### Running RAiSD
list=$(cat "$vcf/scaff_name.txt")

for i in $list
do
RAiSD -n $i -I $vcf/$i.recode.vcf -D -P -O -A 0.99 -s -t -R 
done











.
