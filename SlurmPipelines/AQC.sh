#!/bin/bash
FILE=$1
OUTNAME=${FILE/R1.fastq.gz/.pbs}
echo "#!/bin/bash
#SBATCH -p quick
#SBATCH --job-name=AQC${FILE/.bam/}
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-2:00:00



after.py -1 $FILE -2 ${FILE/R1.fastq/R2.fastq} -g ./good/ -b ./bad/ -r ./afterQC/

">AQC$OUTNAME

sbatch AQC$OUTNAME
mv AQC$OUTNAME batchScripts/.