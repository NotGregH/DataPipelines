

import subprocess
from csv import writer
import argparse as argp

def capture_file_names(folder):
	file_list = subprocess.check_output(['ls',folder])
	file_list = file_list.split("\n")
	filesPeaks = []
	for i in file_list:
		if i != '':
			if ".bam" in i and '.bai' not in i:
				if i not in filesPeaks:
					filesPeaks.append(i)
	return filesPeaks 

def generateCoverageScripts(folder,bedFile,name,outFolder):
	files = capture_file_names(folder)

	print("Files Captured", files)

	for i in files:
		pbsFile = i+'_DepthAt_'+name+'.pbs'
		with open(pbsFile,'w') as o:
			write = writer(o,delimiter=" ",lineterminator='\n')
			
			write.writerow(('#!/bin/bash',''))
			write.writerow(('#$','-q','all.q'))
			write.writerow(('#$','-l','h_vmem=4G'))
			write.writerow(('#$','-cwd'))
			write.writerow(('#$','-j', 'y'))
			write.writerow(('#$','-R','y'))
			write.writerow(('#$','-N',name+'_Coverage_'+i))

			write.writerow(("module","purge"))
			write.writerow(('module','load','samtools/1.5/gcc.4.4.7'))
			write.writerow(("module","load","python/2.7.8/gcc.4.4.7"))
			write.writerow(("mkdir",outFolder))

			write.writerow(("python","GetCoverage.py",
				"-b",i,
				"-l",bedFile,
				"-n",name,
				"-f",folder,
				"-o",outFolder
				))
		o.close()

		subprocess.call(['qsub',pbsFile])

	return


if __name__ == '__main__':
	parser = argp.ArgumentParser()
	parser.add_argument('-f','--folder')
	parser.add_argument('-b','--bedFile')
	parser.add_argument('-n','--name')
	parser.add_argument('-o','--outFolder')
	arg = parser.parse_args()
	f = arg.folder
	name = arg.name
	b = arg.bedFile
	o = arg.outFolder
	generateCoverageScripts(f,b,name,o)
