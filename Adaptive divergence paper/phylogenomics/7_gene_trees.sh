############################################################
######################## Gene trees ########################
############################################################
#!/bin/bash
source activate iqtree
mkdir -p 6gene_trees

ls 4align_mafft/*.fa | xargs -P 4 -I {} bash -c '
  f="{}"
  base=$(basename "$f" .fa)
  iqtree2 -s "$f" -m MFP -bb 1000 -T 20 -pre "6gene_trees/$base"
'

