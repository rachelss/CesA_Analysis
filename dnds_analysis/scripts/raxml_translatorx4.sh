#!/bin/bash
#SBATCH --job-name="RAxML"
#SBATCH --time=96:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node 
#SBATCH --mail-user="kerrina.whelan@uri.edu" #CHANGE to your email
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

module purge

module load RAxML/8.2.12-intel-2019b-hybrid-avx2

ALIGNMENT=/data/schwartzlab/kawhelan/data/CesA_transx_alignment_aliview_edit.fasta #Put the path to your alignment file here
OUTPUT=CesA_transx_alignment_aliview_edit.output #Give your output file a name relevant to the alignment
P=36 #Change to the number of processors (20 or 36)

raxmlHPC -s $ALIGNMENT -n $OUTPUT -m GTRGAMMA -T $P -f a -p $RANDOM -N 100 -x $RANDOM

