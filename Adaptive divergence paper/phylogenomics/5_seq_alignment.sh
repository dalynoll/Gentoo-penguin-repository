############################################################################
######################## Align sequences with MAFFT ########################
############################################################################

#!/bin/bash
source activate mafft
mkdir -p 4align_mafft

ls 3single_copy/*.fa | \
xargs -P 30 -I {} bash -c '
  f="{}"
  b=$(basename "$f" .fa)
  echo "Alineando $f ..."
  mafft --auto "$f" > "4align_mafft/${b}_aligned.fa"
'

