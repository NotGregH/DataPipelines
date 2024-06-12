"""
Date - 11/14/18


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
	pairedFiles = []
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
				i = i.replace('.fastq','')
				i = i.replace('.fq','')
				i = i.replace('.gz','')
				if i not in filesPeaks:
					filesPeaks.append(i)
				else:
					if i not in pairedFiles:
						pairedFiles.append(i)

	return filesPeaks , gzip, pairedFiles


def generatePreProcessScripts(folder,genomeDir,p,q,m,intron,uniqueMap):
	'''
	This function creates the PBS 
	scripts and submits them to the que
	chosen by the arguments. 


	'''

	FASTQ_Files, gzip, pairedFiles = capture_file_names(folder+'/fastq/')

	print("Fastqs Captured", FASTQ_Files)
	print("Paired End Files : ", pairedFiles)
	singleFiles = []
	for i in FASTQ_Files:
		if i not in pairedFiles:
			singleFiles.append(i)
		
	print("Single End Files : ", singleFiles)
	for i in FASTQ_Files:

		pbsFile = i+'_Process.pbs'
		with open(pbsFile,'w') as o:
			write = writer(o,delimiter=" ",lineterminator='\n')
			
			# Setting pbs script variables
			write.writerow(('#!/bin/bash',''))
			write.writerow(('#$','-q',q))
			write.writerow(('#$','-l','h_vmem='+m))
			write.writerow(('#$','-cwd'))
			write.writerow(('#$','-j', 'y'))
			write.writerow(('#$','-R','y'))
			write.writerow(('#$','-N','Seq_'+i+'_Processing'))
			# Loading all required modules
			write.writerow(('module','load','trim_galore/0.4.1'))
			write.writerow(('module','load','FastQC/0.11.4/java.1.8.0_20'))
			write.writerow(('module','load','STAR/2.6.1b/gcc.4.9.2'))
			write.writerow(('module','load','picard-tools/2.3.0/java.1.8.0_20'))
			write.writerow(('module','load','samtools/1.5/gcc.4.4.7'))
			write.writerow(('module','load','bedtools2/2.26.0/gcc.4.4.7'))
			# Setting up directory
			write.writerow(("cd", folder))
			write.writerow(("mkdir", "fastqTrim"))
			write.writerow(("mkdir", "fastqTrim/fastqc"))
			write.writerow(("mkdir","BAM"))
			write.writerow(("mkdir","BAM/alignmentLogs"))
			write.writerow(("mkdir","BAM/DupStatFiles"))
			write.writerow(("mkdir","BAM/unmarkedBams"))
			write.writerow(("mkdir","BAM/samtoolsStats"))

			
			if i in pairedFiles:
				# Trim
				write.writerow(("cd","fastq"))
				write.writerow(("gunzip",i+"_R1.fastq.gz"))
				write.writerow(("gunzip",i+"_R2.fastq.gz"))
				write.writerow(("mv",i+"_R1.fastq",i+"_1.fastq"))
				write.writerow(("mv",i+"_R2.fastq",i+"_2.fastq"))
				write.writerow(("cd","./../"))			
				write.writerow(("trim_galore", "--paired","-o", "fastqTrim",
					"./fastq/"+i+"_1.fastq", "./fastq/"+i+"_2.fastq" ))
				write.writerow(("cd", "./fastqTrim/"))
				#fastqc Trim
				write.writerow(("mv",i+"_1_val_1.fq", i+"_1.fq"))
				write.writerow(("mv",i+"_2_val_2.fq", i+"_2.fq"))
				write.writerow(("mkdir", "fastqc"))
				write.writerow(("fastqc","-t", p, "-o", "fastqc", i+"_1.fq"))
				write.writerow(("fastqc","-t", p, "-o", "fastqc", i+"_2.fq"))
				write.writerow(("cd", "./../"))
				#align
				if uniqueMap == 'Y':
					write.writerow(('STAR','--alignIntronMax',intron,'--outFilterMultimapNmax',
						'1','--outFileNamePrefix','./BAM/'+i,'--runThreadN', p,'--outSAMtype',
						'BAM', 'SortedByCoordinate','--genomeDir',genomeDir,'--readFilesIn',
						'./fastqTrim/'+i+'_1.fq','./fastqTrim/'+i+'_2.fq'))
				if uniqueMap == 'N':
					write.writerow(('STAR','--alignIntronMax',intron,
						'--outFileNamePrefix','./BAM/'+i,'--runThreadN', p,'--outSAMtype',
						'BAM', 'SortedByCoordinate','--genomeDir',genomeDir,'--readFilesIn',
						'./fastqTrim/'+i+'_1.fq','./fastqTrim/'+i+'_2.fq'))			

			# Single Read Proccessing
			if i in singleFiles:
				# Trim
				write.writerow(("trim_galore", "-o", "fastqTrim", "./fastq/"+i+".fastq" ))
				write.writerow(("cd", "./fastqTrim/"))
				#fastqc Trim
				write.writerow(("mv",i+"_trimmed.fq", i+".fq"))
				write.writerow(("fastqc","-t", p, "-o", "fastqc", i+".fq"))
				write.writerow(("cd", "./../"))
				#align
				if uniqueMap == 'Y': 
					write.writerow(('STAR','--alignIntronMax',intron,'--outFilterMultimapNmax',
						'1','--outFileNamePrefix','./BAM/'+i,'--runThreadN', p,'--outSAMtype', 
						'BAM', 'SortedByCoordinate','--genomeDir',genomeDir,
						'--readFilesIn','./fastqTrim/'+i+'.fq'))			
				if uniqueMap == 'N': 
					write.writerow(('STAR','--alignIntronMax',intron,
						'--outFileNamePrefix','./BAM/'+i,'--runThreadN', p,'--outSAMtype', 
						'BAM', 'SortedByCoordinate','--genomeDir',genomeDir,
						'--readFilesIn','./fastqTrim/'+i+'.fq'))		
			# Bam Process
			write.writerow(('cd','BAM'))
			write.writerow(('java', '-jar', '$(which', 'picard.jar)','MarkDuplicates', 
				'M='+i+'_dupstats.txt',"REMOVE_DUPLICATES=true", 'I='+i+'Aligned.sortedByCoord.out.bam',
				 'O='+i+'_NoDups.bam'))
			write.writerow(('samtools', 'index', '-@', p, i+'_NoDups.bam', i+'_NoDups.bam.bai'))
			write.writerow(('samtools','idxstats',i+'_NoDups.bam','>',i+'_ST.txt'))

			#Cleaning Directory
			write.writerow(('mv',i+'*Log.*','alignmentLogs/.'))
			write.writerow(('mv',i+'*SJ.out.tab','alignmentLogs/.'))
			write.writerow(('mv',i+'*.sortedByCoord*','unmarkedBams/.'))
			write.writerow(('mv',i+'*dupstats*','DupStatFiles/.'))
			write.writerow(('mv',i+'*_ST.txt','samtoolsStats/.'))


		o.close()
		subprocess.call(['qsub',pbsFile])

	return 


if __name__ == '__main__':
	parser = argp.ArgumentParser()
	parser.add_argument('-f','--folder')
	parser.add_argument('-g','--genomeDir')
	parser.add_argument('-p','--processors',default='16')
	parser.add_argument('-m','--mem',default='60G')
	parser.add_argument('-q','--que',default='all.q')
	parser.add_argument('-i','--intronMaxSize',default='1')
	parser.add_argument('-u','--uniqueMap',default='Y')
	arg = parser.parse_args()
	folder = arg.folder
	genomeDir = arg.genomeDir	
	p = arg.processors
	q = arg.que
	m = arg.mem
	i = arg.intronMaxSize
	u = arg.uniqueMap
	generatePreProcessScripts(folder,genomeDir,p,q,m,i,u)
