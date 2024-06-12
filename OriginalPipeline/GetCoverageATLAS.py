

import subprocess
import argparse as argp
from csv import writer
import os


def getCoverage(bam):
 	total = 0
	with open("../m1Kd_Run1/macroh2a1_1_ISOR.bed",'r') as f:
		with open("../m1Kd_Run1/m1_Coverage","w") as o:
			write = writer(o,delimiter="\t")
			for line in f:
				line = line.strip().split("\t")
				loc = line[0]+":"+line[1]+"-"+line[2]
				command = "samtools view -F 0x40 ./../m1Kd_Run1/alignedTrim_Unique/"+bam+"_Marked_NOMT.bam " +loc + " | cut -f 6 | grep -v 'N' | wc -l"
				coverage = subprocess.check_output(command,shell=True)
				line.append(int(coverage))
				write.writerow(line)
				total = total + int(coverage)


			o.close()
		f.close()
		
	print(bam)
	print(total)
	return





if __name__ == '__main__':
	parser = argp.ArgumentParser()
	parser.add_argument('-b','--bam')
	arg = parser.parse_args()
	bam = arg.bam
	getCoverage(bam)
