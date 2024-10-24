#!/bin/bash
FILE=$1

echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=trimgalore${FILE/fastq.gz/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./../slurmout/R-%x.%j.out
#SBATCH --time=0-04:00:00
mkdir ./TrimedFastq
trim_galore --nextseq 10 --stringency 5 -o ./TrimedFastq -j 10 --paired --gzip $FILE ${FILE/R1.fastq.gz/R2.fastq.gz}
">trimgalore$OUTNAME

sbatch trimgalore$OUTNAME
mv trimgalore$OUTNAME ./batchScripts/.

