# ATACseq_bulk
QC and pipeline of ATACseq bulk of FOUNDIN

### Typical pipeline below:

```
#!/bin/bash
FILENAME=$1
# start like this:
# sbatch --cpus-per-task=20 --mem=100g --mail-type=BEGIN,END --time=8:00:00 mapping.sh SAMPLEID (e.g. PPMI3658_d0)

echo "$FILENAME"_R1_001.fastq.gz

# Check read quality, output is .html

module load fastqc/0.11.9

fastqc /data/CARD/FOUNDINtemp/ATACseq/fastq/"$FILENAME"_R1_001.fastq.gz

fastqc /data/CARD/FOUNDINtemp/ATACseq/fastq/"$FILENAME"_R2_001.fastq.gz

# Map reads to GRCh38 human genome with bowtie2, final output is .bam 

module load bowtie/2-2.4.1

bowtie2 -p 10 --local -x /data/reedx/FBn_ATACseq/GRCh38/Sequence/Bowtie2Index/genome -1 /data/CARD/FOUNDINtemp/ATACseq/fastq/"$FILENAME"_R1_001.fastq.gz -2 /data/CARD/FOUNDINtemp/ATACseq/fastq/"$FILENAME"_R2_001.fastq.gz -S /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam

echo 'wc -l aligned'
wc -l /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam

echo 'wc -l chrM'
grep chrM /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam | wc -l

echo 'wc -l chrUn'
grep chrUn /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam | wc -l

sed '/chrM/d;/random/d' < /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam > /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_removedchrs.sam

module load samtools/1.3.1

samtools view -bu /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_removedchrs.sam | samtools sort -o /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_sorted.bam

rm /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_aligned.sam

samtools index /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_sorted.bam

module load samtools/0.1.19

samtools rmdup /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_sorted.bam /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.bam
 
samtools view /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.bam -o /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.sam
	
samtools index /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.bam

echo 'wc -l filtered'
wc -l /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.sam

rm /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_filtered.sam

rm /data/CARD/FOUNDINtemp/ATACseq/mapped_reads/"$FILENAME"_sorted.bam

# Compile sequencing statistics from slurm*.out files into an Excel file
```
