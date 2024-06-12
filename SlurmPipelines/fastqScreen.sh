#!/bin/bash
FILE=$1
OUTDIR=$2
OUTNAME=${FILE/"_R1.fastq.gz"/.pbs}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=fastqScreen${FILE/"_R1.fastq.gz"/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-10:00:00
fastq_screen --conf /gs/gsfs0/home/ghamilto/Greg/2022_ChIPData/rawFastqs/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf

mkdir $OUTDIR
gunzip $FILE
gunzip ${FILE/R1.fastq/R2.fastq}
fastq_screen --threads 16 --outdir $OUTDIR ${FILE/.gz/} ${FILE/R1.fastq.gz/R2.fastq}

">fastqScreen$OUTNAME

sbatch fastqScreen$OUTNAME
mv fastqScreen$OUTNAME batchScripts/.

