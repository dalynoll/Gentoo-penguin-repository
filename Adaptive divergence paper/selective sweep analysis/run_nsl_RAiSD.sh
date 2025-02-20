list=$(cat "scaff_name.txt")

for i in $list;
do vcftools  --vcf  $VCF_FILE  --chr $i  --recode --recode-INFO-all --out  $i&
done
wait```



###### Running RAiSD

```bash
list=$(cat "$vcf/scaff_name.txt")

for i in $list
do
RAiSD -n ${i} -I $vcf/${i}.vcf -D -P -O -A 0.99 -s -t -R 
done



###### Running nSL
```bash
list=$(cat "$vcf/scaff_name.txt")

for i in $list
