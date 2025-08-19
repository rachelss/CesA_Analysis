#!/bin/bash
#SBATCH -c 1  # Request 1 CPU (codeml is single-threaded)
#SBATCH --mem=15G  
#SBATCH -p cpu  # Use the 'cpu' partition
#SBATCH --qos=long  # Use long queue (enables >48h wall time)
#SBATCH -t 48:00:00  
#SBATCH -o /work/pi_aroberts_uri_edu/slurm/sitemodels-%j.out  # Write SLURM output to log file (%j = job ID)

# Move into working directory

cd /work/pi_aroberts_uri_edu/cesa_dnds/paml/paml4.9j/src/  

# Define model numbers
# NSsites = 0 (one-ratio), 1 (neutral), 2 (positive selection), 7 (beta), 8 (beta+omega)

MODELS=(0 1 2 7 8)

# Input files
 
ALIGNMENT="CesA_transx_alignment_aliview_edit2.phy"
TREE="RAxML_bestTree.CesA_transx_alignment_aliview_edit.output"
TEMPLATE="codeml-sites.ctl"

# Output directories

OUTDIR="site_model_results"
mkdir -p $OUTDIR slurm

# Run each model
# Loop over each model

for model in "${MODELS[@]}"; do
    echo "Running codeml with NSsites = $model"

    # Define custom control file and output path

    CTL="ctl_site_model${model}.ctl"
    OUT="site_model${model}.txt"

    # Copy template and replace key lines

    cp $TEMPLATE $CTL
    sed -i "s|seqfile = .*|seqfile = $ALIGNMENT|" $CTL
    sed -i "s|treefile = .*|treefile = $TREE|" $CTL
    sed -i "s|outfile = .*|outfile = $OUT|" $CTL
    sed -i "s|NSsites = .*|NSsites = $model|" $CTL

    # Run codeml with model-specific control file

    ./codeml $CTL

done

echo "All site models completed."

