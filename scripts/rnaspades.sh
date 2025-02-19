#!/bin/bash
#SBATCH -J SPAdes
#SBATCH -t 100:00:00
#SBATCH -n 36
#SBATCH -N 1
#SBATCH -p schwartzlab
#SBATCH --mem=100g
#SBATCH --export=NONE
#SBATCH -o SPAdes.out
#SBATCH -e SPAdes.err

module purge
module load SAMtools/1.14-GCC-11.2.0
module load Bowtie2/2.4.4-GCC-11.2.0
module load seqtk/1.3-GCC-11.2.0
module load SPAdes/3.15.3-GCC-11.2.0
module load GCCcore/11.2.0

export PATH=$PATH:/data/schwartzlab/cbreusing/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/schwartzlab/cbreusing/lib

# Change file and lib names for other species and SRA accessions
file=Fontinalis_antipyretica
lib1=ERR2040966

bowtie2 --very-sensitive-local -p ${SLURM_CPUS_ON_NODE} -x ../CESA_references.fasta -1 ../../filtered_reads/${file}_${lib1}_R1_clean.fastq.gz -2 ../../filtered_reads/${file}_${lib1}_R2_clean.fastq.gz | samtools view -bS -h -@ ${SLURM_CPUS_ON_NODE} - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.sorted.bam
samtools index ${file}.sorted.bam

samtools view -h -@ ${SLURM_CPUS_ON_NODE} -F4 ${file}.sorted.bam > ${file}.CESA.sam
cut -f1 ${file}.CESA.sam | sort | uniq > ${file}.CESA_ids.lst
seqtk subseq ../../filtered_reads/${file}_${lib1}_R1_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib1}_R1.fastq
seqtk subseq ../../filtered_reads/${file}_${lib1}_R2_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib1}_R2.fastq
rm ${file}.CESA.sam
rnaspades.py -k 21,33,55,77 -t ${SLURM_CPUS_ON_NODE} --pe1-1 ${file}_CESA_${lib1}_R1.fastq --pe1-2 ${file}_CESA_${lib1}_R2.fastq -o ${file}_CESA
perl ../RenameContigs.pl ${file}_CESA/transcripts.fasta ${file}_CESA ${file}_CESA.fasta
rm -r ${file}_CESA
