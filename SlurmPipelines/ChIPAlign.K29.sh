#!/bin/bash
FILE=$1
OUTDIR=$2
GENOME=$3
BLACKLIST=$4
OUTNAME=${FILE/_R1.trim.fastq.gz/.pbs}
BAM=${FILE/_R1.trim.fastq.gz/.Align.bam}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=BWA${FILE/.fastq.gz/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-05:00:00

cd ../
mkdir $OUTDIR
cd $OUTDIR
mkdir STFS
mkdir STIS
mkdir Picard
bwa mem -M -k 29 -t 40 $GENOME ./../rawFastqs/$FILE ./../rawFastqs/${FILE/R1.trim/R2.trim} | samtools view -b -@ 40 -h -S - | samtools sort -@ 40 -o $BAM -
samtools index $BAM
samtools flagstat $BAM  > ./STFS/${BAM/.bam/.FS.txt}
samtools idxstat $BAM > ./STIS/${BAM/.bam/.IS.txt}
picard CollectInsertSizeMetrics I=$BAM O=./Picard/${BAM/.bam/.INS.txt}  H=./Picard/${BAM/.bam/.INS.pdf} W=1000
picard MarkDuplicates M=./Picard/${BAM/.bam/.dup.txt} REMOVE_DUPLICATES=false I=$BAM O=${BAM/Align/Marked}
samtools index ${BAM/Align/Marked}
samtools flagstat ${BAM/Align/Marked}  > ./STFS/${BAM/Align.bam/Marked.FS.txt}
samtools idxstat ${BAM/Align/Marked} > ./STIS/${BAM/Align.bam/Marked.IS.txt}
samtools view --verbosity 0 -b -h -@ 30 -L $BLACKLIST -U ${BAM/Align/BLF} ${BAM/Align/Marked}
samtools index ${BAM/Align/BLF}
samtools flagstat ${BAM/Align/BLF}  > ./STFS/${BAM/Align.bam/BLF.FS.txt}
samtools idxstat ${BAM/Align/BLF} > ./STIS/${BAM/Align.bam/BLF.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Align/BLF} O=./Picard/${BAM/Align.bam/BLF.INS.txt}  H=./Picard/${BAM/Align.bam/BLF.INS.pdf} W=1000
samtools view --verbosity 0 -b -h -@ 30 -q 10 -f 2 -F 256 ${BAM/Align/BLF} > ${BAM/Align/Final}
samtools index ${BAM/Align/Final}
samtools flagstat ${BAM/Align/Final}  > ./STFS/${BAM/Align.bam/Final.FS.txt}
samtools idxstat ${BAM/Align/Final} > ./STIS/${BAM/Align.bam/Final.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Align/Final} O=./Picard/${BAM/.bam/.Final.INS.txt}  H=./Picard/${BAM/.bam/.Final.INS.pdf} W=1000
mkdir bsCov
mkdir bsCovBLF
mkdir bsCovFin
BAMscale scale -t 20 -f -o bsCov -n ${BAM/.bam/} --bam $BAM
BAMscale scale -t 20 -f -o bsCovBLF -n ${BAM/Align.bam/BLF} --bam ${BAM/Align/BLF}
BAMscale scale -t 20 -f -o bsCovFin -n ${BAM/Align.bam/Final} --bam ${BAM/Align/Final}
">BWA.Align.$OUTNAME

sbatch BWA.Align.$OUTNAME
mv BWA.Align.$OUTNAME batchScripts/.

