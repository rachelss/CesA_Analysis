#!/bin/bash

#SBATCH -c 1

#SBATCH --mem=15G

#SBATCH -p cpu

#SBATCH --qos=long

#SBATCH -t 60:00:00

#SBATCH -o /work/pi_aroberts_uri_edu/slurm/branchloop-%j.out



# Set working directory

cd /work/pi_aroberts_uri_edu/cesa_dnds/paml/paml4.9j/src/


echo "✅ SLURM script started on: $(date)" > slurm_debug.log

echo "Current directory: $(pwd)" >> slurm_debug.log

echo "Files here:" >> slurm_debug.log

ls -l >> slurm_debug.log



# Run test codeml call

## ./codeml codeml-branch.ctl >> slurm_debug.log 2>&1  # Test if needed



echo "✅ codeml finished on: $(date)" >> slurm_debug.log

# Input files

TREE="RAxML_bestTree.CesA_transx_alignment_aliview_edit.output"

ALIGNMENT="CesA_transx_alignment_aliview_edit2.phy"

CTL_TEMPLATE="codeml-branch.ctl"

OUTPUT_DIR="branch_results"



mkdir -p $OUTPUT_DIR



# Extract unique taxon names from tree (roughly parsed)

LABELS=$(grep -oP "\)\K[^):]+" $TREE | sort | uniq)



# Loop through each label and run codeml

for LABEL in $LABELS; do

    echo "🧪 Running codeml for branch: $LABEL"



    # Make tree with current branch labeled as foreground

    TREE_MOD="tree_${LABEL}.nwk"

    sed "s/\(${LABEL}\)/\1#1/" $TREE > $TREE_MOD



    # Create ctl file for this branch

    CTL_MOD="ctl_${LABEL}.ctl"

    cp $CTL_TEMPLATE $CTL_MOD



    sed -i "s|seqfile = .*|seqfile = $ALIGNMENT|" $CTL_MOD

    sed -i "s|treefile = .*|treefile = $TREE_MOD|" $CTL_MOD

    sed -i "s|outfile = .*|outfile = $OUTPUT_DIR/out_${LABEL}.txt|" $CTL_MOD



    # Run codeml

    ./codeml $CTL_MOD

done



echo "✅ All branch models completed."

