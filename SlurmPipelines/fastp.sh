#!/bin/bash
FILE=$1
OUTNAME=${FILE/R1.fastq.gz/.pbs}
echo "#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=fastp${FILE/.fastq.qz/}
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --output=./slurmout/R-%x.%j.out
#SBATCH --time=0-06:00:00

mkdir fastp

fastp -w 16 -i $FILE -I ${FILE/R1.fastq/R2.fastq} --dedup --dont_overwrite -o ${FILE/R1.fastq/R1.fp.fastq} -O ${FILE/R1.fastq/R2.fp.fastq} -c -g -j ./fastp/${FILE/.fastq.gz/.json} -h ./fastp/${FILE/.fastq.gz/.html} -R ${FILE/R1.fastq.gz/} --detect_adapter_for_pe --trim_poly_g --trim_poly_x --dup_calc_accuracy 6
fastqc -t 12 -o ./fastfpqc/ ${FILE/R1.fastq/R1.fp.fastq} ${FILE/R1.fastq/R2.fp.fastq} 
">fastp$OUTNAME

sbatch fastp$OUTNAME
mv fastp$OUTNAME batchScripts/.