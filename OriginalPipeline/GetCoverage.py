

import subprocess
import argparse as argp
from csv import writer
import os


def getCoverage(bam,bed,folder,outFolder,name):
 	total = 0
	with open(bed,'r') as f:
		with open(outFolder+name+bam.replace('.bam','.txt'),"w") as o:
			write = writer(o,delimiter="\t")
			for line in f:
				outline = []
				line = line.strip().split("\t")
				loc = line[0]+":"+line[1]+"-"+line[2]
				outline.append(loc)
				command = "samtools view " +folder+ bam + " " + loc + " | cut -f 6 | grep -v 'N' | wc -l"
				coverage = subprocess.check_output(command,shell=True)
				outline.append(int(coverage))
				write.writerow(outline)
				total = total + int(coverage)

		o.close()
	f.close()

	return





if __name__ == '__main__':
	parser = argp.ArgumentParser()
	parser.add_argument('-b','--bam')
	parser.add_argument('-l','--bed')
	parser.add_argument('-f','--folder')
	parser.add_argument('-o','--outFolder')
	parser.add_argument('-n','--name')
	arg = parser.parse_args()
	bam = arg.bam
	bed = arg.bed
	folder = arg.folder
	outFolder = arg.outFolder
	name = arg.name
	getCoverage(bam,bed,folder,outFolder,name)
