#!/bin/bash
FILE=$1
OUTDIR=$2
GENOME=$3
BLACKLIST=$4
OUTNAME=${FILE/_R1.fp.fastq.gz/.pbs}
BAM=${FILE/_R1.fp.fastq.gz/Aligned.sortedByCoord.out.bam}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=STAR${FILE/fastq.gz/}
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
STAR --readFilesCommand gunzip --alignIntronMax 1 --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNoverReadLmax 0.05 --outFileNamePrefix ${FILE/_R1.fp.fastq.gz/} --runThreadN 20 --outSAMtype BAM SortedByCoordinate --genomeDir $GENOME --readFilesIn ./../rawFastqs/$FILE ./../rawFastqs/${FILE/R1.fp/R2.fp}
samtools index $BAM
samtools flagstat $BAM  > ./STFS/${BAM/.bam/.FS.txt}
samtools idxstat $BAM > ./STIS/${BAM/.bam/.IS.txt}
picard CollectInsertSizeMetrics I=$BAM O=./Picard/${BAM/.bam/.INS.txt}  H=./Picard/${BAM/.bam/.INS.pdf} W=1000
picard MarkDuplicates M=./Picard/${BAM/.bam/.dup.txt} REMOVE_DUPLICATES=false I=$BAM O=${BAM/Aligned.sortedByCoord.out/Marked}
samtools index ${BAM/Aligned.sortedByCoord.out/Marked}
samtools flagstat ${BAM/Aligned.sortedByCoord.out/Marked}  > ./STFS/${BAM/Aligned.sortedByCoord.out.bam/Marked.FS.txt}
samtools idxstat ${BAM/Aligned.sortedByCoord.out/Marked} > ./STIS/${BAM/Aligned.sortedByCoord.out.bam/Marked.IS.txt}
samtools view --verbosity 0 -b -h -@ 30 -L $BLACKLIST -U ${BAM/Aligned.sortedByCoord.out/BLF} ${BAM/Aligned.sortedByCoord.out/Marked}
samtools index ${BAM/Aligned.sortedByCoord.out/BLF}
samtools flagstat ${BAM/Aligned.sortedByCoord.out/BLF}  > ./STFS/${BAM/Aligned.sortedByCoord.out/BLF.FS.txt}
samtools idxstat ${BAM/Aligned.sortedByCoord.out/BLF} > ./STIS/${BAM/Aligned.sortedByCoord.out.bam/BLF.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Aligned.sortedByCoord.out/BLF} O=./Picard/${BAM/Aligned.sortedByCoord.out.bam/BLF.INS.txt}  H=./Picard/${BAM/Aligned.sortedByCoord.out.bam/BLF.INS.pdf} W=1000
samtools view --verbosity 0 -b -h -@ 30 -q 10 -f 2 -F 256 ${BAM/Aligned.sortedByCoord.out/BLF} > ${BAM/Aligned.sortedByCoord.out/Final}
samtools index ${BAM/Aligned.sortedByCoord.out/Final}
samtools flagstat ${BAM/Aligned.sortedByCoord.out/Final}  > ./STFS/${BAM/Aligned.sortedByCoord.out.bam/Final.FS.txt}
samtools idxstat ${BAM/Aligned.sortedByCoord.out/Final} > ./STIS/${BAM/Aligned.sortedByCoord.out.bam/Final.IS.txt}
picard CollectInsertSizeMetrics I=${BAM/Aligned.sortedByCoord.out/Final} O=./Picard/${BAM/Aligned.sortedByCoord.out.bam/.Final.INS.txt}  H=./Picard/${BAM/Aligned.sortedByCoord.out.bam/.Final.INS.pdf} W=1000
mkdir bsCov
mkdir bsCovBLF
mkdir bsCovFin
BAMscale scale -t 20 -f -o bsCov -n ${BAM/.bam/} --bam $BAM
BAMscale scale -t 20 -f -o bsCovBLF -n ${BAM/Aligned.sortedByCoord.out.bam/BLF} --bam ${BAM/Aligned.sortedByCoord.out/BLF}
BAMscale scale -t 20 -f -o bsCovFin -n ${BAM/Aligned.sortedByCoord.out.bam/Final} --bam ${BAM/Aligned.sortedByCoord.out/Final}
">STAR.Align.$OUTNAME

sbatch STAR.Align.$OUTNAME
mv STAR.Align.$OUTNAME batchScripts/.

