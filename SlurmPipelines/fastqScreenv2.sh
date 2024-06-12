#!/bin/bash
FILE=$1
OUTNAME=${FILE/"_R1.fastq"/.pbs}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=fastqScreen${FILE/"_R1.fastq"/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-10:00:00


fastq_screen --conf ./FastQ_Screen_Genomes/fastq_screen.conf --threads 16 $FILE ${FILE/R1.fastq/R2.fastq}

">fastqScreen$OUTNAME

sbatch fastqScreen$OUTNAME
mv fastqScreen$OUTNAME batchScripts/.

