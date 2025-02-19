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

# Change file and lib named for other species and SRA accessions
file=Ceratodon_purpureus_R40
lib1=SRR4125803
lib2=SRR4125804

ln -s ../../filtered_reads/${file}_${lib1}_R1_clean.fastq.gz ${file}_${lib1}_R1_clean.fastq.gz
ln -s ../../filtered_reads/${file}_${lib1}_R2_clean.fastq.gz ${file}_${lib1}_R2_clean.fastq.gz
ln -s ../../filtered_reads/${file}_${lib2}_R1_clean.fastq.gz ${file}_${lib2}_R1_clean.fastq.gz
ln -s ../../filtered_reads/${file}_${lib2}_R2_clean.fastq.gz ${file}_${lib2}_R2_clean.fastq.gz

bowtie2 --very-sensitive-local -p ${SLURM_CPUS_ON_NODE} -x ../CESA_references.fasta -1 ${file}_${lib1}_R1_clean.fastq.gz,${file}_${lib2}_R1_clean.fastq.gz -2 ${file}_${lib1}_R2_clean.fastq.gz,${file}_${lib2}_R2_clean.fastq.gz | samtools view -bS -h -@ ${SLURM_CPUS_ON_NODE} - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.sorted.bam
samtools index ${file}.sorted.bam

samtools view -h -@ ${SLURM_CPUS_ON_NODE} -F4 ${file}.sorted.bam > ${file}.CESA.sam
cut -f1 ${file}.CESA.sam | sort | uniq > ${file}.CESA_ids.lst
seqtk subseq ${file}_${lib1}_R1_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib1}_R1.fastq
seqtk subseq ${file}_${lib1}_R2_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib1}_R2.fastq
seqtk subseq ${file}_${lib2}_R1_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib2}_R1.fastq
seqtk subseq ${file}_${lib2}_R2_clean.fastq.gz ${file}.CESA_ids.lst > ${file}_CESA_${lib2}_R2.fastq
rm ${file}.CESA.sam
spades.py -k 21,33,55,77 -t ${SLURM_CPUS_ON_NODE} --pe1-1 ${file}_CESA_${lib1}_R1.fastq --pe2-1 ${file}_CESA_${lib2}_R1.fastq --pe1-2 ${file}_CESA_${lib1}_R2.fastq --pe2-2 ${file}_CESA_${lib2}_R2.fastq -o ${file}_CESA
perl ../RenameContigs.pl ${file}_CESA/scaffolds.fasta ${file}_CESA ${file}_CESA.fasta
rm -r ${file}_CESA
