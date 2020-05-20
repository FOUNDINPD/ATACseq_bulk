# ATACseq_bulk
QC and pipeline of ATACseq bulk of FOUNDIN-PD

Xylena Reed May 2020

### Typical pipeline below:

```
#!/bin/bash
FILENAME=$1
WORKINGDIR = /PATH/TO/WORKING/DIR/
# start like this:
# sbatch --cpus-per-task=20 --mem=100g --mail-type=BEGIN,END --time=8:00:00 ATACseq_bulk.sh SAMPLEID

echo "$FILENAME"_R1_001.fastq.gz

# Check read quality, output is .html

module load fastqc/0.11.9

fastqc /$WORKINGDIR/fastq/"$FILENAME"_R1_001.fastq.gz

fastqc /$WORKINGDIR/fastq/"$FILENAME"_R2_001.fastq.gz

# Map reads to GRCh38 human genome with bowtie2, final output is .bam 

module load bowtie/2-2.4.1

bowtie2 -p 10 --local -x /$WORKINGDIR/GRCh38/Sequence/Bowtie2Index/genome -1 /$WORKINGDIR/fastq/"$FILENAME"_R1_001.fastq.gz -2 /$WORKINGDIR/fastq/"$FILENAME"_R2_001.fastq.gz -S /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam

echo 'wc -l aligned'
wc -l /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam

echo 'wc -l chrM'
grep chrM /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam | wc -l

echo 'wc -l random'
grep random /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam | wc -lt

sed '/chrM/d;/random/d;/chrUn/d' < /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam > /$WORKINGDIR/mapped_reads/"$FILENAME"_removedchrs.sam

module load samtools/1.3.1

samtools view -bu /$WORKINGDIR/mapped_reads/"$FILENAME"_removedchrs.sam | samtools sort -o /$WORKINGDIR/mapped_reads/"$FILENAME"_sorted.bam

rm /$WORKINGDIR/mapped_reads/"$FILENAME"_aligned.sam

samtools index /$WORKINGDIR/mapped_reads/"$FILENAME"_sorted.bam

module load samtools/0.1.19

samtools rmdup /$WORKINGDIR/mapped_reads/"$FILENAME"_sorted.bam /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.bam
 
samtools view /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.bam -o /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.sam
	
samtools index /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.bam

echo 'wc -l filtered'
wc -l /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.sam

rm /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.sam

rm /$WORKINGDIR/mapped_reads/"$FILENAME"_sorted.bam

# Compile sequencing statistics from slurm*.out files into an Excel file

module load macs/2.1.2

# Peaks representing Nucleosome Free Regions (NFRs)

macs2 callpeak -t /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.bam -n /$WORKINGDIR/peaks/"$FILENAME"_NFRs --nomodel --shift -100 --extsize 200 --format BAM -g hs

# regular peak calling 

macs2 callpeak --nomodel -B -f BAMPE -t /$WORKINGDIR/mapped_reads/"$FILENAME"_filtered.bam -n /$WORKINGDIR/peaks/"$FILENAME"_orig --nolambda --keep-dup all --gsize hs

```
