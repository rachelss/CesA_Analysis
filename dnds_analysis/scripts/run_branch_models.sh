#!/bin/bash

#SBATCH -c 1

#SBATCH --mem=15G

#SBATCH -p cpu

#SBATCH --qos=long

#SBATCH -t 75:00:00

#SBATCH -o /work/pi_aroberts_uri_edu/slurm/branchloop-%j.out



cd /work/pi_aroberts_uri_edu/cesa_dnds/paml/paml4.9j/src/



ALIGNMENT="CesA_transx_alignment_aliview_edit2.phy"

CTL_TEMPLATE="codeml-branch.ctl"

TREE_DIR="branch_foreground_trees"

OUTPUT_DIR="branch_results"



mkdir -p $OUTPUT_DIR



for TREE_FILE in $TREE_DIR/tree_*.nwk; do

    # Extract label from filename

    LABEL=$(basename "$TREE_FILE" .nwk | sed 's/^tree_//')



    echo "🔄 Running codeml for: $LABEL"



    CTL_COPY="ctl_${LABEL}.ctl"

    cp $CTL_TEMPLATE $CTL_COPY



    sed -i "s|seqfile = .*|seqfile = $ALIGNMENT|" $CTL_COPY

    sed -i "s|treefile = .*|treefile = $TREE_FILE|" $CTL_COPY

    sed -i "s|outfile = .*|outfile = $OUTPUT_DIR/out_${LABEL}.txt|" $CTL_COPY



    ./codeml $CTL_COPY

done



echo "✅ All codeml runs complete."


