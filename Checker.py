#A script to systematically check for each gene whether the SuperTranscript builder worked ok

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.pyplot import cm
import seaborn as sns
from matplotlib import gridspec

##################################################
###### Visualise blocks in SuperTranscript #######
##################################################

def Checker(genome):

	print("Finding list of genes")
	genes=[]
	f = open(sys.argv[1],'r')
	for line in f:
		if('>' in line):
			genes.append((line.split('>')[1]).split("\n")[0])

	metrics = {}
	for gene in genes:
		#print(gene)
		mapping,fraction = FindMetrics(gene)
		metrics[gene] = [mapping,fraction]

	#Fraction of genes where we get a one to one mapping
	mapp_frac = 0
	frac_covered = []
	for key in metrics:
		mapp_frac += metrics[key][0]	
		for key2 in metrics[key][1]:
			frac_covered.append(metrics[key][1][key2])

	mapp_frac = mapp_frac / len(genes)
	frac_covered = np.asarray(frac_covered)


	#Now lets pring some info        
	print("Fraction of genes with a one to one mapping:")
	print(mapp_frac)

	#Plot the distribution
	n, bins, patches = plt.hist(frac_covered,20,normed=1,facecolor='green',alpha=0.75)
	plt.xlabel("Fraction of transcript covered in SuperTranscript")
	plt.ylabel("Frequency")
	plt.title(r'$\mathrm{Histogram\ of\ transcript\ buidling\ of\ ST}$')
	plt.show()
		

def FindMetrics(gene_name):

	#EXTRACT GENE FROM SUPER
	
	#Find gene in genome
	f= open("SuperDuper.fasta","r")
	for line in f:
		if(gene_name in line):
			gene_string=next(f)
			break
	f.close()

	fs= open("Super.fasta","w")
	fs.write(">" + gene_name + "\n")
	fs.write(gene_string)
	fs.close()


	#Match transcripts to super transcript
	print("Producing match to super transcript")
	BLAT_command = "./blat Super.fasta %s supercomp.psl" %(gene_name)
	os.system(BLAT_command)

	#First read in BLAT output:
	Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']
	vData = pd.read_table('supercomp.psl',sep='\t',header = None,names=Header_names,skiprows=5)

	#Make list of transcripts
	transcripts = np.unique(list(vData['Q name']))

	#Metric 1 - If each transcript has a one to one mapping with ST then we should have as many lines as transcripts
	mapping = 0
	if(len(transcripts) == len(vData)): mapping =1

	#Metric 2 - The fraction fo each transcript mapped (i.e the number of bases)
	fraction = {} # Dictonary of fraction mapped
	for i in range(0,len(vData)):		
		if vData.iloc[i,9] in fraction:  fraction[vData.iloc[i,9]] += (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10]) #If key already in dictionary sum fractions
		fraction[vData.iloc[i,9]] = (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10])
		
	return(mapping, fraction)

if __name__=='__main__':
	if(len(sys.argv) != 2):
		print("Checker function requires one argument")
		print("The genome whose super transcripts you wish to quality check")
	else:
		Checker(sys.argv[1])
