#!/bin/bash
#SBATCH -J ReadFilter
#SBATCH -t 24:00:00
#SBATCH -n 18
#SBATCH -N 1
#SBATCH -p schwartzlab
#SBATCH --mem=100g
#SBATCH -o ReadFilter_%a.out
#SBATCH -e ReadFilter_%a.err
#SBATCH --array=1-187%10

module purge
module load FastQC/0.11.9-Java-11 
module load fastp/0.23.2-GCC-11.2.0
module load Trimmomatic/0.39-Java-11
module load SAMtools/1.14-GCC-11.2.0
module load Bowtie2/2.4.4-GCC-11.2.0
module load seqtk/1.3-GCC-11.2.0
module load pigz/2.6-GCCcore-11.2.0

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" files.list)

fastqc -t ${SLURM_CPUS_ON_NODE} ../raw_reads/${file}_1.fastq.gz
fastqc -t ${SLURM_CPUS_ON_NODE} ../raw_reads/${file}_2.fastq.gz

fastp --in1 ../raw_reads/${file}_1.fastq.gz --out1 ${file}_gtrim_1.fq --in2 ../raw_reads/${file}_2.fastq.gz --out2 ${file}_gtrim_2.fq -A -Q -L -g -y -w ${SLURM_CPUS_ON_NODE}
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${SLURM_CPUS_ON_NODE} -phred33 ${file}_gtrim_1.fq ${file}_gtrim_2.fq ${file}_paired_R1.fq ${file}_singles_R1.fq ${file}_paired_R2.fq ${file}_singles_R2.fq ILLUMINACLIP:/data/schwartzlab/cbreusing/adaptors/IlluminaBGISEQ.fa:2:30:10:2:True SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75
bowtie2 -p ${SLURM_CPUS_ON_NODE} -x /data/schwartzlab/cbreusing/databases/contaminants -1 ${file}_paired_R1.fq -2 ${file}_paired_R2.fq | samtools view -bS -h -@ ${SLURM_CPUS_ON_NODE} - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.bowtie2.cont.sorted.bam
samtools view -@ ${SLURM_CPUS_ON_NODE} -f12 ${file}.bowtie2.cont.sorted.bam > ${file}.cont.unmapped.sam
cut -f1 ${file}.cont.unmapped.sam | sort | uniq > ${file}.cont.unmapped_ids.lst
seqtk subseq ${file}_paired_R1.fq ${file}.cont.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_paired_R2.fq ${file}.cont.unmapped_ids.lst > ${file}_R2_clean.fastq
pigz -p ${SLURM_CPUS_ON_NODE} ${file}_R1_clean.fastq
pigz -p ${SLURM_CPUS_ON_NODE} ${file}_R2_clean.fastq

rm ${file}_gtrim_*.fq
rm ${file}_paired_R*.fq
rm ${file}_singles_R*.fq
rm ${file}.cont.unmapped.sam
