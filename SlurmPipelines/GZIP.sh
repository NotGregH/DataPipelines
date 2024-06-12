#!/bin/bash
FILE=$1
OUTNAME=${FILE/.fastq/.pbs}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=gzip${FILE/fastq.gz/}
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-05:00:00

gzip $FILE
">gzip.$OUTNAME

sbatch gzip.$OUTNAME
mv gzip.$OUTNAME batchScripts/.

