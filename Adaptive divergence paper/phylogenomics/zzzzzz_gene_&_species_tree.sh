############################################################
######################## Gene trees ########################
############################################################
#!/bin/bash
source activate iqtree
mkdir -p 5_gene_trees_iqtree

ls 4_align_mafft_2269/*.fa | xargs -P 4 -I {} bash -c '
  f="{}"
  base=$(basename "$f" .fa)
  iqtree2 -s "$f" -m MFP -bb 1000 -T 8 -pre "5_gene_trees_iqtree/$base"
'

cat *.treefile > all_gene_trees.tre

###############################################################
######################## Species trees ########################
###############################################################

#!/bin/bash
mkdir -p 6_SpeciesTree_astral

java -jar /data2/daly/gentoo/filogenias/cds/programas/ASTRAL/Astral/astral.5.7.8.jar \
  -i list_gene_trees.txt \
  -o species_tree.tre
