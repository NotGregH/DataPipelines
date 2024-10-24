#!/bin/bash
FILE=$1
OUTNAME=${FILE/_val_1.fq.gz/.pbs}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=fastqc${FILE/_val_1.fq.gz/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --output=./../slurmout/R-%x.%j.out
#SBATCH --time=0-02:00:00
mkdir ./fastqc

mv $FILE ${FILE/_val_1.fq/.trim.fastq}
mv ${FILE/R1_val_1/R2_val_1} ${FILE/R1_val_1.fq/R2.trim.fastq}
fastqc -t 12 -o ./fastqc/ ${FILE/_val_1.fq/.trim.fastq} ${FILE/R1_val_1.fq/R2.trim.fastq}

">fastqc.$OUTNAME

sbatch fastqc.$OUTNAME
mv fastqc.$OUTNAME ./../batchScripts/.

