#!/bin/bash
#SBATCH -J raxml
#SBATCH -t 100:00:00
#SBATCH -n 36
#SBATCH -N 1
#SBATCH -p uri-cpu
#SBATCH --mem=100g
#SBATCH -o raxml.out
#SBATCH -e raxml.err

eval "$(conda shell.bash hook)"
conda activate raxml
# Run this command for the full CESA alignment and the reference CESA alignment
raxmlHPC-PTHREADS -s CESA_exon.aln -n CESA_exon.tre -m GTRGAMMA -T ${SLURM_CPUS_ON_NODE} -f a -p 123 -N autoMRE -x 123
conda deactivate