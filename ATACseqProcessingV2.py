
"""
Date - 10/10/19


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


"""

__author__ =  'Gregory A. Hamilton'
__version__ = '0.2.0'
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


def generatePreProcessScripts(folder,genomeDir,p,q,m,uniqueMap):
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
			write.writerow(('#$','-N','ATAC_'+i+'_Processing'))
			# Loading all required modules
			write.writerow(('module','load','trim_galore/0.4.1'))
			write.writerow(('module','load','FastQC/0.11.4/java.1.8.0_20'))
			write.writerow(('module','load','STAR/2.6.1b/gcc.4.9.2'))
			write.writerow(('module','load','picard-tools/2.3.0/java.1.8.0_20'))
			write.writerow(('module','load','samtools/1.5/gcc.4.4.7'))
			write.writerow(('module','load','deeptools/3.1.0/python.2.7.8'))
			write.writerow(('module','load','MACS2/2.1.0-update/python.2.7.8'))
			write.writerow(('module','load','genrich/0.5/gcc.4.4.7'))
			write.writerow(('module','load','bedtools2/2.26.0/gcc.4.4.7'))
			# Setting up directory
			write.writerow(("cd", folder))
			write.writerow(("mkdir","BAM"))
			write.writerow(("mkdir","fastqTrim"))
			write.writerow(("mkdir","BAM/alignmentLogs"))
			write.writerow(("mkdir","BAM/DupStatFiles"))
			write.writerow(("mkdir","BAM/unmarkedBams"))
			write.writerow(("mkdir","BAM/samtoolsStats"))
			write.writerow(("mkdir","fastqTrim"))
			write.writerow(("mkdir","fragDist"))

			
			if i in pairedFiles:
				if uniqueMap == 'Y':
					write.writerow(("trim_galore", "--paired","-o", "fastqTrim",
						"./fastq/"+i+"_R1.fastq.gz", "./fastq/"+i+"_R2.fastq.gz" ))
					write.writerow(("cd", "./fastqTrim/"))
					#fastqc Trim
					write.writerow(("mv",i+"_R1_val_1.fq.gz", i+"_1.fq.gz"))
					write.writerow(("mv",i+"_R2_val_2.fq.gz", i+"_2.fq.gz"))
					write.writerow(("mkdir", "fastqc"))
					write.writerow(("fastqc","-t", p, "-o", "fastqc", i+"_1.fq.gz"))
					write.writerow(("fastqc","-t", p, "-o", "fastqc", i+"_2.fq.gz"))
					write.writerow(("cd", "./../"))
					write.writerow(('STAR','--readFilesCommand','gunzip','-c','--alignIntronMax','1','--outFilterMultimapNmax',
						'1','--outFileNamePrefix','./BAM/'+i,'--runThreadN', p,'--outSAMtype',
						'BAM', 'SortedByCoordinate','--genomeDir',genomeDir,'--readFilesIn',
						'./fastqTrim/'+i+'_1.fq.gz','./fastqTrim/'+i+'_2.fq.gz','--outFilterScoreMinOverLread', '0.4' ,
						'--outFilterMatchNminOverLread', '0.4'))		
			# Bam Process
			write.writerow(('cd','BAM'))
			write.writerow(('java', '-jar', '$(which', 'picard.jar)','MarkDuplicates','M='+i+'_dupstats.txt',"REMOVE_DUPLICATES=true", 'I='+i+'Aligned.sortedByCoord.out.bam','O='+i+'_NoDups.bam'))
			write.writerow(('samtools', 'index', '-@', p, i+'_NoDups.bam', i+'_NoDups.bam.bai'))
			write.writerow(('samtools','idxstats',i+'_NoDups.bam','>',i+'_ST.txt'))

			#Cleaning Directory
			write.writerow(('mv',i+'*Log.*','alignmentLogs/.'))
			write.writerow(('mv',i+'*SJ.out.tab','alignmentLogs/.'))
			write.writerow(('mv',i+'*.sortedByCoord*','unmarkedBams/.'))
			write.writerow(('mv',i+'*dupstats*','DupStatFiles/.'))
			write.writerow(('mv',i+'*_ST.txt','samtoolsStats/.'))

			# Filtering BAM File
			write.writerow(('samtools','view','-@',p,'-b','-F','4',i+'_NoDups.bam','1','2','3','4','5','6','7',
				'8','9','10','11','12' ,'13' ,'14', '15' ,'16' ,'17' ,'18' ,'19' ,'20' ,'21','22','>', i+'_Final.bam'))
			write.writerow(('samtools', 'index', '-@', p, i+'_Final.bam', i+'_Final.bam.bai'))

			write.writerow(("mkdir","PreQCbam"))
			write.writerow(("mkdir","DeepToolsLogs"))
			write.writerow(('mv',i+'*_NoDups_*','PreQCbam/.'))
			# Calculating Frag Distrobution 
			write.writerow(('samtools','view','-@','16',i+'_Final.bam','|','awk','\'$9>0\'','|','cut','-f','9','|','sort','|','uniq','-c','|','sort','-b','-k2,2n','|','sed','-e',"s/^[ \\t]*//",'>','../FragDist/'+i+'.txt'))
			write.writerow(("mkdir","../MACS2"))
			write.writerow(("mkdir","../BigWigs"))
			## Creating BigWigs
			write.writerow(('bamCoverage','-b',i+'_Final.bam','-o','../BigWigs/'+i+'.bw','-bs=1','--normalizeUsing=RPGC','-p='+p,'--exactScaling'))
			write.writerow(('bamCoverage','-b',i+'_Final.bam','-o','../BigWigs/'+i+'_50.bw','-bs=50','--normalizeUsing=RPGC','-p='+p,'--exactScaling'))
			write.writerow(('bamCoverage','-b',i+'_Final.bam','-o','../BigWigs/'+i+'_100.bw','-bs=100','--normalizeUsing=RPGC','-p='+p,'--exactScaling'))
			## Calling Peaks
			write.writerow(('macs2','callpeak','--nomodel','--extsize','100','-t',i+'_Final.bam','-g','hs','--keep-dup','all','-n','../MACS2/'+i,'--shift','-50','--call-summits'))
			write.writerow(('macs2','callpeak','-t',i+'_Final.bam','-g','hs','--keep-dup','all','-n','../MACS2/'+i+"_BAMPE","-f","BAMPE"))
			write.writerow(("mkdir","../Genrich"))
			write.writerow(("samtools", "sort", "-@",p, "-n" ,i+"_Final.bam", ">", i+"_Final.sort.bam"))
			write.writerow(("Genrich","-a","100","-f","../Genrich/"+i+".log","-t",i+"_Final.sort.bam","-j","-y","-v","-o","../Genrich/"+i+"_100.narrowPeak","-E","../PeakFiles/ENCFF001TDO.bed"))
			write.writerow(("Genrich","-a","200","-P","-f","../Genrich/"+i+".log","-j","-y","-v","-o","../Genrich/"+i+"_200.narrowPeak","-E","../PeakFiles/ENCFF001TDO.bed"))
			write.writerow(("Genrich","-P","-f","../Genrich/"+i+".log","-t","-j","-y","-v","-o","../Genrich/"+i+"_20.narrowPeak","-E","../PeakFiles/ENCFF001TDO.bed"))
			write.writerow(('cd','unmarkedBams'))
			write.writerow(('samtools', 'index', '-@',p,i+'Aligned.sortedByCoord.out.bam', i+'Aligned.sortedByCoord.out.bam.bai'))
			write.writerow(('samtools','idxstats',i+'Aligned.sortedByCoord.out.bam','>',i+'_ST.txt'))
			write.writerow(('cd','../'))
			write.writerow(('mkdir','ATAC_QC_Metrics'))
			write.writerow(('echo','\'Reads In Genrich Peaks:\'','>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('bedtools', 'sort', '-i', "../Genrich/"+i+"_20.narrowPeak",'|', 'bedtools', 'merge', 
				'-i', 'stdin', '|', 'bedtools', 'intersect', '-u', '-nonamecheck', '-a', i+'_Final.bam', '-b', 'stdin', '-ubam',
				'|', 'samtools', 'view', '-c','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('echo','\\n','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('echo','Reads In Genrich 100 Peaks:','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('bedtools', 'sort', '-i', "../Genrich/"+i+"_100.narrowPeak",'|', 'bedtools', 'merge', 
				'-i', 'stdin', '|', 'bedtools', 'intersect', '-u', '-nonamecheck', '-a', i+'_Final.bam', '-b', 'stdin', '-ubam',
				'|', 'samtools', 'view', '-c','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))

			write.writerow(('echo','\\n','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('echo','Reads In Genrich 200 Peaks:','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('bedtools', 'sort', '-i', "../Genrich/"+i+"_200.narrowPeak",'|', 'bedtools', 'merge', 
				'-i', 'stdin', '|', 'bedtools', 'intersect', '-u', '-nonamecheck', '-a', i+'_Final.bam', '-b', 'stdin', '-ubam',
				'|', 'samtools', 'view', '-c','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))

			write.writerow(('echo','\\n','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('echo','Reads In MACS2 Peaks:','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('bedtools', 'sort', '-i', '../MACS2/'+i+'_peaks.narrowPeak','|', 'bedtools', 'merge', 
				'-i', 'stdin', '|', 'bedtools', 'intersect', '-u', '-nonamecheck', '-a', i+'_Final.bam', '-b', 'stdin', '-ubam',
				'|', 'samtools', 'view', '-c','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))

			write.writerow(('echo','\\n','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('echo','Total Reads:','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))
			write.writerow(('samtools','view','-c',i+'_Final.bam','>>','ATAC_QC_Metrics/'+i+'_ReadsInPeaks.txt'))

			write.writerow(('bedtools', 'bamtobed', '-bedpe', '-i', i+'_Final.sort.bam', '|', 'awk', 'BEGIN{OFS=\\t}{print $1,$2,$4,$6,$9,$10}','|', 
				'grep', '-v', 'chrM', '|', 'sort', '|', 'uniq', '-c', '|', 'awk', 
				'\'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf' ,'\"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n', 'mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}', '>','ATAC_QC_Metrics/'+i+".pbc.qc"))

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
	parser.add_argument('-u','--uniqueMap',default='Y')
	arg = parser.parse_args()
	folder = arg.folder
	genomeDir = arg.genomeDir	
	p = arg.processors
	q = arg.que
	m = arg.mem
	u = arg.uniqueMap
	generatePreProcessScripts(folder,genomeDir,p,q,m,u)
