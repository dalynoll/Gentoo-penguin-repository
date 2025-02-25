#!/bin/bash

source activate agat
anotacion=/data2/daly/gentoo/filogenias/gff_curados/gff
genomas=/data2/daly/gentoo/fasta/gentooref
cds=/data2/daly/gentoo/filogenias/gff_curados/1fasta_cds_curated

file="13A2
15A2
18A1
19A1
22A1
24A1
CRO05
CRO07
CRO08
G4
G5
G6
G7
GEN04
GEN05
GEN08
GEN09
GEN10
GSA10
GSA12
GSA14
GSA15
GSA7
GSA8
GSA9
GTP06
GTP07
GTP12
MQ9209
MQ9312
P0048
P0051
P0052
P0056
P0058
P0065
P0167
P0169
P0170
P0173
P0196
P0198
P0200
P0402
P0404
P0405
P0407
P0408
P0409
P0412
P1884
P1885
P1887
P1889
P1893
P1896
P2134
P2136
P2143
PAM10
PAM1
PAM3
PAM5
PAM7
adel
chins"

for i in $file
do
agat_sp_extract_sequences.pl -g $anotacion/gentoo_genomic_filtered.gff -f $genomas/${i}_maskedPygo.fa -o $cds/${i}_cds.fa -t cds
done

#############################################################################################
# for South Georgia 
#############################################################################################
#!/bin/bash

source activate agat
anotacion=/data2/daly/gentoo/filogenias/gff_curados/gff
genomas=/data2/daly/gentoo/fasta/sgeorgia
cds=/data2/daly/gentoo/filogenias/gff_curados/1fasta_cds_curated

agat_sp_extract_sequences.pl -g $anotacion/southgeorgias_genomic_filtered.gff -f $genomas/southgeorgia_assembly.fa -o $cds/sgeorgia_cds.fa -t cds -p


#############################################################################################
# simplify_headers.sh
#############################################################################################
#!/usr/bin/env bash
#
# Script: simplify_headers.sh
# Use: ./simplify_headers.sh /path/to/input_dir /path/to/output_dir
#
#
set -e

if [ "$#" -ne 2 ]; then
    echo "Uso: $0 reheader"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Create outpot directory
mkdir -p "$OUTPUT_DIR"

# Loop for each fasta in the Input dir
for fasta in "$INPUT_DIR"/*.fa; do
    BASENAME=$(basename "$fasta")
    echo "Processing: $BASENAME"
    sed -r 's/^>.*gene=(gene-[^ ]+).*/>\1/' "$fasta" > "$OUTPUT_DIR/$BASENAME"
done


