#!/bin/bash
FILE=$1
OUTDIR=$2
OUTNAME=${FILE/_R1.fp.fastq.gz/.pbs}
BAM=${FILE/_R1.fp.fastq.gz/.Align.bam}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=BWA${FILE/_R1.fp.fastq.gz/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-10:00:00

cd ../
mkdir $OUTDIR
cd $OUTDIR
mkdir STFS
mkdir STIS
mkdir Picard
bwa-mem2 mem -M -t 40 ./../../Genomes/hg19mm10v3/bwamm10v2/mm10 ./../rawFastqs/$FILE ./../rawFastqs/${FILE/R1.fp/R2.fp} | samtools view -@ 30 -b -h -S -F 4 - | samtools sort -@ 10 -o $BAM -
samtools index $BAM
samtools flagstat $BAM  > ./STFS/${BAM/.bam/.FS.txt}
samtools idxstat $BAM > ./STIS/${BAM/.bam/.IS.txt}
picard MarkDuplicates M=./Picard/${BAM/.bam/.dup.txt} REMOVE_DUPLICATES=false I=$BAM O=${BAM/Align/Marked}
samtools index ${BAM/Align/Marked}
samtools flagstat ${BAM/Align/Marked}  > ./STFS/${BAM/Align.bam/Marked.FS.txt}
samtools idxstat ${BAM/Align/Marked} > ./STIS/${BAM/Align.bam/Marked.IS.txt}
samtools view -b -h -@ 30 -L ./../Blacklist.bed -U ${BAM/Align/BLF} ${BAM/Align/Marked}
samtools index ${BAM/Align/BLF}
samtools flagstat ${BAM/Align/BLF}  > ./STFS/${BAM/Align.bam/BLF.FS.txt}
samtools idxstat ${BAM/Align/BLF} > ./STIS/${BAM/Align.bam/BLF.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Align/BLF} O=./Picard/${BAM/.bam/.INS.txt}  H=./Picard/${BAM/.bam/.INS.pdf} W=1000
samtools view -b -h -@ 30 -q 55 -f 2 -F 256 ${BAM/Align/BLF} > ${BAM/Align/Final}
samtools index ${BAM/Align/Final}
samtools flagstat ${BAM/Align/Final}  > ./STFS/${BAM/Align.bam/Final.FS.txt}
samtools idxstat ${BAM/Align/Final} > ./STIS/${BAM/Align.bam/Final.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Align/Final} O=./Picard/${BAM/.bam/.Final.INS.txt}  H=./Picard/${BAM/.bam/.Final.INS.pdf} W=1000
">BWA.Align.$OUTNAME

sbatch BWA.Align.$OUTNAME
mv BWA.Align.$OUTNAME batchScripts/.

