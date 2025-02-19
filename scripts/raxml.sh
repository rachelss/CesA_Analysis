#!/bin/bash
#SBATCH -J raxml
#SBATCH -t 100:00:00
#SBATCH -n 36
#SBATCH -N 1
#SBATCH -p schwartzlab
#SBATCH --export=NONE
#SBATCH --mem=100g
#SBATCH -o raxml.out
#SBATCH -e raxml.err

module purge
module load RAxML/8.2.12-intel-2019b-hybrid-avx2

# Run this command for the full CESA alignment and the reference CESA alignment
raxmlHPC -s CESA_exon.aln -n CESA_exon.tre -m GTRGAMMA -T ${SLURM_CPUS_ON_NODE} -f a -p 123 -N 100 -x 123
