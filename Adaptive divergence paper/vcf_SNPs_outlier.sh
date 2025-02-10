

cat *_shared_windows.txt |awk '{print $1"\t"$4"\t"$5}' |grep -v '^chr' >gentoo_all_outlier_windows.bed

bedtools sort -i gentoo_all_outlier_windows.bed > gentoo_all_outlier_windows_sorted.bed

bedtools merge -i gentoo_all_outlier_windows_sorted.bed > gentoo_all_outlier_windows_sorted_merged.bed



vcftools --gzvcf $vcf/gentoo.vcf.gz --keep gentoo62_toVCFoutlier.txt --bed gentoo_all_outlier_windows_sorted_merged.bed --recode --out gentoo_outliers
 ### After filtering, kept 31234 out of a possible 9158384 Sites


