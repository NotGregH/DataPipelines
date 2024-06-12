"""
Date - 5/25/18


Writen for the Gamble Lab @ 
Albert Einstein College of Medicine


This script is intended to be used to
generate pbs scripts to process all next gen
sequencing data.  

# Pre Alignement
The script will first capture the names 
of all fastq files to be processed. The files
will then be trimmed and place the trimmed
fastqs in /fastqTrim/ .  The trimmed fastqs
will then be analyzed with fastqc. 

# Alignment
The fastqs will then be aligned with STAR aligner. 
The output bam files will be sorted by coordinate and
only contain uniquely aligned reads. 

# Bam Proccessing
The duplicates in the bam files will then be flagged
using picard tools. The marked bam files will then be indexed. 

# Arguments

-f - Path/to/the/base/project/folder
-g - Path/to/the/star/genome/
-pE - Y or N  (Y is for paired end reads) (Default = Y)
-p - # of processors ( default = 16 )
-q - The que the jobs are to be submitted to ( default = all.q)
-m - The amount of memmory requested (defualt = 60G)


# Plans for future updates:
1) Enable RNA sequencing alignment (right now the setting do not allow 
for alignment accross splice sites)

"""

__author__ =  'Gregory A. Hamilton'
__version__ = '0.0.2'
__email__ = 'ghamilto@mail.einstein.yu.edu'

###################
# Required Modules
###################

import sys
import argparse as argp
from csv import writer
import subprocess


def capture_file_names(folder):
	'''
	Goes to a folder and captures all unique
	fastq file names. For paired end it only
	captures the name once.
	'''
	file_list = subprocess.check_output(['ls',folder])
	file_list = file_list.split("\n")
	filesPeaks = []
	gzip = 'N'
	for i in file_list:
		if i != '':
			if ".fastq" in i or '.fq' in i:
				i = i.replace('_R1.fastq','')
				i = i.replace('_R2.fastq','')
				i = i.replace('_R1.fq','')
				i = i.replace('_R2.fq','')
				i = i.replace('_1.fastq','')
				i = i.replace('_2.fastq','')
				i = i.replace('_1.fq','')
				i = i.replace('_2.fq','')
				i = i.replace('.fastq','')
				i = i.replace('.fq','')
				i = i.replace('.gz','')
				if i not in filesPeaks:
					filesPeaks.append(i)
	return filesPeaks , gzip


def generatePreProcessScripts(folder,genomeDir,pE,p,q,m,intron):
	'''
	This function creates the PBS 
	scripts and submits them to the que
	chosen by the arguments. 


	'''

	FASTQ_Files, gzip = capture_file_names(folder+'/fastq/')

	print("Fastqs Captured", FASTQ_Files)


	for i in FASTQ_Files:

		pbsFile = i+'_ProcessBowtie2.pbs'
		with open(pbsFile,'w') as o:
			write = writer(o,delimiter=" ",lineterminator='\n')
			
			# Setting pbs script variables
			write.writerow(('#!/bin/bash',''))
			write.writerow(('#$','-q',q))
			write.writerow(('#$','-l','h_vmem='+m))
			write.writerow(('#$','-cwd'))
			write.writerow(('#$','-j', 'y'))
			write.writerow(('#$','-R','y'))
			write.writerow(('#$','-N','Seq_'+i+'_ProcessingBowtie2'))
			# Loading all required modules
			write.writerow(('module','load','trim_galore/0.4.1'))
			write.writerow(('module','load','FastQC/0.11.4/java.1.8.0_20'))
			write.writerow(('module','load','bowtie2/2.3.3.1'))
			write.writerow(('module','load','picard-tools/2.3.0/java.1.8.0_20'))
			write.writerow(('module','load','samtools/1.5/gcc.4.4.7'))
			# Setting up directory
			write.writerow(("cd", folder))
			write.writerow(("mkdir","Bowtie2alignment"))
			write.writerow(("mkdir","Bowtie2alignment/DupStatFiles"))
			write.writerow(("mkdir","Bowtie2alignment/unmarkedBams"))
			write.writerow(("mkdir","Bowtie2alignment/samtoolsStats"))
			write.writerow(("mkdir","Bowtie2alignment/"+i))



			write.writerow(("cd", "./Bowtie2alignment/"+i))
			write.writerow(('bowtie2','-X2000','-p',p,'-x',genomeDir,'-1','./../../fastqTrim/'+i+'_1.fq','-2','./../../fastqTrim/'+i+'_2.fq','-S','../'+i+'.sam','>>','./../'+i+'_bowtieLog.txt'))
			write.writerow(("cd","../"))
			write.writerow(("samtools","view","-S","-b",i+".sam",">",i+".unsorted.bam"))
			write.writerow(("samtools","sort",i+".unsorted.bam","-o",i+".bam"))	

			# Bam Process
			write.writerow(('java', '-jar', '$(which', 'picard.jar)','MarkDuplicates', 
				'M='+i+'_dupstats.txt',"REMOVE_DUPLICATES=true", 'I='+i+'.bam',
				 'O='+i+'_NoDups.bam'))
			write.writerow(('samtools', 'index', '-@', p, i+'_NoDups.bam', i+'_NoDups.bam.bai'))
			write.writerow(('samtools', '-h', i+'_NoDups.bam','|',''))
			write.writerow(('samtools','idxstats',i+'_NoDups.bam','>',i+'NoDups_ST.txt'))
			write.writerow(('samtools','idxstats',i+'.bam','>',i+'Dups_ST.txt'))

			#Cleaning Directory
			write.writerow(('mv',i+'.bam','unmarkedBams/.'))
			write.writerow(('mv',i+'*dupstats*','DupStatFiles/.'))
			write.writerow(('mv',i+'*_ST.txt','samtoolsStats/.'))

		o.close()
		subprocess.call(['qsub',pbsFile])

	return 


if __name__ == '__main__':
	parser = argp.ArgumentParser()
	parser.add_argument('-f','--folder')
	parser.add_argument('-g','--genomeDir')
	parser.add_argument('-p','--processors',default='20')
	parser.add_argument('-pE','--pairedEnd',default='Y')
	parser.add_argument('-m','--mem',default='60G')
	parser.add_argument('-q','--que',default='all.q')
	parser.add_argument('-i','--intronMaxSize',default='1')
	arg = parser.parse_args()
	folder = arg.folder
	genomeDir = arg.genomeDir	
	pE = arg.pairedEnd
	p = arg.processors
	q = arg.que
	m = arg.mem
	i = arg.intronMaxSize
	generatePreProcessScripts(folder,genomeDir,pE,p,q,m,i)
