############################################################################
######################## Align sequences with MAFFT ########################
############################################################################

#!/bin/bash
source activate mafft
mkdir -p 4_align_mafft

ls 3_seq_genes/*.fa | \
xargs -P 30 -I {} bash -c '
  f="{}"
  b=$(basename "$f" .fa)
  echo "Alineando $f ..."
  mafft --auto "$f" > "4_align_mafft/${b}_aligned.fa"
'

