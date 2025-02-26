########################################
########### BWA reads mapping ##########
########################################

#!/bin/bash

WD=/path/to/output/
FQ=/path/to/cleanedFQ/
REF=/path/to/GCA_010090195.1_BGI_Ppap.V1_genomic.fna

~/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 10 -M -R '@RG\tID:sample\tLB:lib1\tPL:ILLUMINA\tPU:unit1\tSM:sample' \
$REF \
$FQ/sample_filtered_1P.fq \
$FQ/sample_filtered_2P.fq > $WD/sample_clean.sam

samtools view -S -b sample_clean.sam > sample_clean.bam



####################################
########### BAM filtering ##########
####################################


  # # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # #
# # # # # # properly paired reads and read quality ≥ 10 # # # # # #
  # # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # #

#!/bin/bash
files="sample1
sample2
sample3
..."

for sample in $files;
do
samtools view -bh -q 10 -f 0x2 -@ 8  ${sample}_clean.bam > ${sample}_q10f0x2.bam
done



  # # # # # # # ## # # # # ## # # # # ## # # # # # # # # #
# # # # # # Sort and mark/remove duplitated reads # # # # # #
  # # # # # # # ## # # # # ## # # # # ## # # # # # # # # #
samtools sort -o ${sample}_sorted.bam \
-@ 20 ${sample}_q10f0x2.bam

picard MarkDuplicates I=${sample}_sorted.bam O=${sample}_sorted_dedup.bam \
METRICS_FILE=${sample}_sorted_dedup_metrics.txt VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true CREATE_MD5_FILE=true TAGGING_POLICY=All ASSUME_SORT_ORDER=coordinate



  # # # # # # ## # # # # ## # # # #
# # # # # # Read realign # # # # # #
  # # # # # # ## # # # # ## # # # #
  
gatk3 -nt 6 \
      -T RealignerTargetCreator \
      -R /path/to/ref/GCA_010090195.1_BGI_Ppap.V1_genomic.fna \
      -I ${sample}_sorted_dedup.bam \
      -o ${sample}_realn.intervals \
      -allowPotentiallyMisencodedQuals -S LENIENT

gatk3 -Xmx64G \
      -T IndelRealigner \
      -R /path/to/ref/GCA_010090195.1_BGI_Ppap.V1_genomic.fna \
      -I ${sample}_sorted_dedup.bam \
      -targetIntervals ${sample}_realn.intervals \
      -o ${sample}_sorted_realign.bam \
      -allowPotentiallyMisencodedQuals -S LENIENT --generate_md5








