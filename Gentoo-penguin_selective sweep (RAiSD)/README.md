
# **Workflow of selective sweep analysis using RAiSD**


## 1. RAiSD
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
