#!/bin/bash
#SBATCH -J SRA
#SBATCH -t 200:00:00
#SBATCH -n 18
#SBATCH -N 1
#SBATCH -p schwartzlab
#SBATCH --mem=100g
#SBATCH -o SRA.out
#SBATCH -e SRA.err

module purge
module load SRA-Toolkit/3.0.3-gompi-2022a
module load pigz/2.7-GCCcore-11.3.0

cp /data/schwartzlab/jezell_greenp/downloaded_data/Downloaded_Data/entodon_seductrix/CL100124498_L01_read_1.fastq.gz Entodon_seductrix_CL100124498_1.fastq.gz
cp /data/schwartzlab/jezell_greenp/downloaded_data/Downloaded_Data/entodon_seductrix/CL100124498_L01_read_2.fastq.gz Entodon_seductrix_CL100124498_2.fastq.gz

cp /data/schwartzlab/jezell_greenp/downloaded_data/Downloaded_Data/hypnum_curvifolium/CL100092064_L01_read_1.fastq.gz Hypnum_curvifolium_CL100092064_1.fastq.gz
cp /data/schwartzlab/jezell_greenp/downloaded_data/Downloaded_Data/hypnum_curvifolium/CL100092064_L01_read_2.fastq.gz Hypnum_curvifolium_CL100092064_2.fastq.gz

# files.list contains a list of species and SRA accessions to use for analysis
IFS=$'\n'

for line in `cat files.list`
do
echo $line > acc.txt
ACC=`cut -f2 acc.txt`
SPEC=`cut -f1 acc.txt`
prefetch --max-size 500G $ACC
fasterq-dump --skip-technical --split-3 --threads ${SLURM_CPUS_ON_NODE} $ACC
rm -r $ACC
pigz -p ${SLURM_CPUS_ON_NODE} ${ACC}_1.fastq
pigz -p ${SLURM_CPUS_ON_NODE} ${ACC}_2.fastq
mv ${ACC}_1.fastq.gz ${SPEC}_${ACC}_1.fastq.gz
mv ${ACC}_2.fastq.gz ${SPEC}_${ACC}_2.fastq.gz
done

