
##############################################################################
# Project Title: Toolkit for piRNA identification in small RNA-seq libraries
#
# Program name: "piRNA_smallRNAseq_analysis.py 
#
# Author: Natalia Llonga
#
# Data: 2018
#
# Python version:  Python 3.6.4 under macOS 10.13.4
#
##############################################################################


'''
Biological problem:

Piwi-interacting RNAs (piRNAs) are small non-coding RNA (sncRNA) of 26 to 31 nucleotides (nt) 
in length. Their system is an evolutionarily conserved mechanism involved in the control of 
transposable elements and maintenance of genomic stability, especially in germ-line cells and 
in early embryo stages. However, relevant particularities, both in mechanism and function, exist 
across species among metazoans.

The identification of piRNAs starts with the novo sequencing small RNA libraries. This means a high 
number of reads that must be analyzed, taking lot of time and spend a huge quantity of RAM memory, 
to complete the analysis. When the analysis must be done in a great number of libraries the problem 
increase, since the study must be done for each library separately. The libraries of RNAseq contains 
a high number of reads that are used to analyze, step by step, usually using shell, python and R.

I propose to create a script in order to facilitate the work, and reduce the time used with each step, 
obtaining as a result, a list of putative piRNAs aligned in the reference genome for each library analyzed. 
I create a script that grouped all of the commands and processes used in this part of the analysis. 
 
In order to identify the piRNAs, I started with libraries of small RNA seq in fastq file format. From fastq files, 
fasta files from each library were obtained and then their reads distribution was visualized. Next, from created 
fasta files those reads that correspond to piRNA fraction (26-31 nucleotides of length) are selected. microRNAs 
(that are previously annotated, and in fasta file), and low complexity sequences were removed. Next reads were 
collapsed in order to obtain unique sequences, that finally were mapped into a reference genome.
The resulting files can be analyzed in R in order to classify the piRNAs by their loci, and also to study their 
expression comparing between different libraries. 

To develop my script I used data from 22 small RNA libraries corresponding to 11 developmental stages from the 
cockroach Blattella germanica (Llonga, et al., under revision)




DEPENDENCES/REQUERIMENTS:

 - 	Bowtie
 -	Bowtie2
 -	Samtools
 -	Fastx_collapser
 -	List of miRNAs in fasta format
 -	Reference genome in fasta format (for example Blattella germanica genome)
 -	All libraries of small RNA-seq to analyse


 '''

##############################################################################
#########################			Libraries			######################

import glob, os, sys 
from Bio import SeqIO
import pandas, pylab

##############################################################################
#########################			Functions			######################

def fastq_to_fasta(file):
	''' This function converts the fastq files to fasta files
	fastq > fasta
	'''
	fileout=str(file)+".fasta"
	SeqIO.convert(file, "fastq", fileout, "fasta")
	print("Fastq to Fasta Done!")


def fasta_select_seqs(file2, min_length, max_length):
	''' This function is used to select reads by their length, from a fasta file
	and returns other fasta file with the sequences 
	fasta > fasta
	'''
	fileout=str(file2)+"26-31.fasta"
	selected_sequences=[]
	for record in SeqIO.parse(file2, "fasta"):
		if len(record.seq)>= min_length and len(record.seq)<= max_length:
			selected_sequences.append(record)
	SeqIO.write(selected_sequences, fileout, "fasta")
	print("Select reads by lenght Done!")

def plot_selected_seqs(file3):
	'''This function is used to show the sequence length distribution 
	from a fasta file
	fasta > jpeg
	'''
	fileout=str(file3)+".jpeg"
	sizes = [len(rec) for rec in SeqIO.parse(file3, "fasta")]
	pylab.hist(sizes, bins=20)
	pylab.title(" \n %i sequences \nLengths %i to %i" \
		% (len(sizes), min(sizes), max(sizes)))
	pylab.xlabel("Sequence length (bp)")
	pylab.ylabel("Count")
	#pylab.show()
	pylab.savefig(fileout)

	print("Histograms done!")

def remove_miRNAs(file3):
	'''
	This function is used to remove miRNA reads from multifasta file, using Bowtie2
	fasta > SAM > BAM > Fastq '''
	COMMAND = ("bowtie2 -L 18 -N 0 -p 2 -x miRNAs_bowtie2 -f "+str(file3)+" -S "+str(file3)+".SAM") #Align reads into a fasta of miRNAs. If the name of fasta file is different, please change "miRNAs_bowtie2" by the name of the file to use.  
	print(COMMAND)
	COMMAND2 = ("samtools view -S "+str(file3)+".SAM -f4 > "+str(file3)+"nomiRNA.SAM") # Select reads that don't align to miRNA sequences
	COMMAND3 = ("samtools view -bS "+str(file3)+"nomiRNA.SAM > "+str(file3)+"nomiRNA.BAM") #SAM to BAM 
	COMMAND4 = ("samtools bam2fq "+str(file3)+"nomiRNA.BAM > "+str(file3)+"nomiRNA.fastq") #BAM to FASTQ
	os.system(COMMAND)
	os.system(COMMAND2)
	os.system(COMMAND3)
	os.system(COMMAND4)
	print ("miRNAs removed")

def remove_low_complexity(file5):
	'''
	Function to remove the fraction of low complexity sequences
	fasta > fasta (with low complexity sequences) and fasta (without low complexity sequences)
	'''
	COMMAND = ("perl TBr2_duster.pl -i "+str(file5)) #TBr2_duster.pl from NGS Toolbox. 
	os.system(COMMAND)
	print("low complexity sequences removed")

def collapse_reads(file6):
	'''
	Function to obtain unique sequences and their counts
	fasta > fasta
	'''
	COMMAND5 = ("fastx_collapser -i "+str(file6)+" -o "+str(file6)+"_collapsed.fasta")
	os.system(COMMAND5)
	print("Collapse sequences done")


def genome_align(file7, genome):
	'''
	Function to align fasta file into a reference genome. Parameters: all best alignments only one mismatch.
	fasta > SAM > BAM
	'''
	COMMAND = ("bowtie -f -a --best --strata -v 1 /Genome/"+genome+" -f "+file7+" --sam "+file7+".SAM" )
	COMMAND2 = ("samtools view -bS "+file7+".SAM > "+file7+".BAM")
	os.system(COMMAND)
	os.system(COMMAND2)
	print("Alignment into a reference genome performed")

##############################################################################
#########################			Program			##########################	

#0 Indexed genome
command_index_bowtie = input("You need an indexed genome. Do you want create one now with Bowtie?(y/n): ")
if command_index_bowtie == "y":
	fa = input("Please, insert the name of the fasta file: ")
	if fa == "Bgermanica.scaffolds.fa":
		os.system("bowtie-build /Genome/Bgermanica.scaffolds.fa /Genome/Bgermanica_bowtie")
	else: 
		os.system("bowtie-build "+fa+" "+fa+"_bowtie")
if command_index_bowtie =="n":
	print("Perfect!")
#0 Indexed miRNA
command2_index_bowtie = input("Do you want create other indexed file with Bowtie2 to miRNA.fasta ?(y/n): ")
if command2_index_bowtie == "y":
	fa2 = input("Please, insert the name of the fasta file: ")
	if fa2 == "miRNAs_gffread.fasta":
		os.system("bowtie2-build miRNAs_gffread.fasta miRNAs_bowtie2")
	else: 
		os.system("bowtie-build "+fa2+" "+fa2+"_bowtie")
if command2_index_bowtie =="n":
	print("Perfect!")


#1 The first step is obtain fasta files from fastq
print("Starting the process ... ")
#os.chdir("/Users/nataliallonga/Desktop/FinalProjectBioinf") #Establish the working directory

for file in glob.glob("*.fastq"):
	print (file)
	fastq_to_fasta(file)

#2 Plot some graphs in order to  check the information

for file2 in glob.glob("*fastq.fasta"):
	print(file2)
	print(plot_selected_seqs(file2))

#3 Select reads of interest. I use piRNA fraction that corresponds between 26 to 31 nucleotides of length

for file2 in glob.glob("*fastq.fasta"):
	fasta_select_seqs(file2, 26, 31)

#4 Plot some graphs in order to  check the information

for file3 in glob.glob("*26-31.fasta"):
	print(file3)
	print(plot_selected_seqs(file3))

#5 Remove miRNA annotated from selected reads. Using Bowtie2

for file3 in glob.glob("*26-31.fasta"):
	print(file3)
	remove_miRNAs(file3)

#6 Obtain fasta file from fastq
for file4 in glob.glob("*nomiRNA.fastq"):
	fastq_to_fasta(file4)

#7 Remove low complexity sequences 

for file5 in glob.glob("*nomiRNA.fastq.fasta"):
	remove_low_complexity(file5)

#7  Collapse reads in order to obtain unique sequences and counts, by using fastx_toolkit

for file6 in glob.glob("*nomiRNA.fastq.fasta.no-dust"):
	collapse_reads(file6)

#8  Align to indexed genome

for file7 in glob.glob("*_collapsed.fasta"):
	Bgermanica = "Bgermanica_bowtie" #I use by default B.germanica genome. If you want to use other, please change the name by other. 
	genome_align(file7, Bgermanica)

###### 
print("Now, the results are OK to study the putative piRNAs ")


