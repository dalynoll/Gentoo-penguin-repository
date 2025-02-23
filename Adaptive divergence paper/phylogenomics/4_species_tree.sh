#####################################################
################## Preliminar tree ##################
#####################################################

## Use catfasta2phyml.pl (https://github.com/nylander/catfasta2phyml) to concatenate all CDS
perl catfasta2phyml.pl *.fa --intersect > cds_10102_gentoo.fa 2> partitions.txt

## Generate a preliminar tree using IQ-Tree

#!/bin/bash
source activate iqtree

iqtree -s cds_10294_gentoo.phy -bb 1000 -T 20 -o adel,chins -m TEST

###### PRELIMINAR TREE OUTPUT:
