### Create plink file including 62 gentoo penguin samples and one outgroup (Chinstrap penguin)
plink --vcf $vcf_gentoo_chin/gentoo62_chin_nomiss.vcf.gz --recode --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### create bfile
plink --file vcf_gentoo_chin --make-bed --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### calculate allele frequencies within sampling localities
plink --bfile vcf_gentoo_chin --freq --missing --within vcf_prun_sinmq_chins.clust --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### compress the frequeny file
gzip vcf_gentoo_chin.frq.strat

### convert plink frequency file to treemix input file (https://github.com/thomnelson/tools/blob/master/plink2treemix.py)
python plink2treemix.py vcf_gentoo_chin.frq.strat.gz vcf_gentoo_chin.treemix.frq.gz

### Run treemix
for m in {1..5}
do
  for i in {1..10}
do
# Generate random seed
s=$RANDOM
echo "Random seed = ${s}"
  echo "Ejecutando TreeMix con m=${m}"
  treemix \
    -i vcf_gentoo_chin.treemix.frq.gz \
    -o run_bootstrap/vcf_gentoo_chin.${i}.${m} \
    -clust vcf_prun_sinmq.clust \
    -m ${m} \
    -global \
    -bootstrap \
    -k 500 \
    -root Chinstrap
    -seed ${s}
  done
done


### Treemix plotting was performing using OptM (Fitak et al 2021)

