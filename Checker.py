#A script to systematically check for each gene whether the SuperTranscript builder worked ok

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.pyplot import cm
import seaborn as sns
from matplotlib import gridspec
import pickle
import time

################################################################
###### Visualise blocks and metricise in SuperTranscript #######
################################################################

def Checker(genome):
	start_time = time.time()
	print("Finding list of genes")
	genes=[]
	f = open(genome,'r')
	for line in f:
		if('>' in line):
			genes.append((line.split('>')[1]).split("\n")[0])
	metrics = {}
	for gene in genes:
		mapping,fraction,anno,compact = FindMetrics(gene)
		metrics[gene] = [mapping,fraction,anno,compact]

	#Fraction of genes where we get a one to one mapping
	mapp_frac = 0
	frac_covered = []
	compactify = []
	for key in metrics:
		mapp_frac += metrics[key][0]
		compactify.append(metrics[key][3])	
		for key2 in metrics[key][1]:
			frac_covered.append(metrics[key][1][key2])

	mapp_frac = mapp_frac / len(genes)
	frac_covered = np.asarray(frac_covered)


	#Now lets pring some info        
	print("Fraction of genes with a one to one mapping:")
	print(mapp_frac)

	#Plot the distribution for metric 2
	n, bins, patches = plt.hist(frac_covered,20,normed=1,facecolor='green',alpha=0.75)
	plt.xlabel("Fraction of transcript covered in SuperTranscript")
	plt.ylabel("Frequency")
	plt.title(r'$\mathrm{Histogram\ of\ transcript\ buidling\ of\ ST}$')
	plt.savefig('Frac_covered.pdf')
	#plt.show()

	#Plot Distribution for Metric 3
	compactified = np.asarray(compactify)
	print(compactified)
	n,bins,patches = plt.hist(compactified,10,facecolor='blue',alpha=0.75)
	plt.xlabel("Sum of bases in Transcripts/ Super Transcript Length")
	plt.ylabel("Frequency")
	plt.savefig('Compactify.pdf')
	#plt.show()
	
	#Save the dicts as pickle file
	pickle.dump(frac_covered,open("frac_covered.pkl","wb"))
	pickle.dump(mapp_frac,open("mapp_frac.pkl","wb"))
	pickle.dump(compactified,open("compactified.pkl","wb"))


	#Now lets save the annotation to file
	fg = open("SuperDuperTrans.gff","w")
	for key in metrics:
		fg.write(metrics[key][2])
	fg.close()
	print("---- %s seconds ----" %(time.time()-start_time))
	
def FindMetrics(gene_name):

	#EXTRACT GENE FROM SUPER
	gene_string=""
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
	BLAT_command = "./blat Super.fasta %s.fasta supercomp.psl" %(gene_name)
	os.system(BLAT_command)

	#First read in BLAT output:
	Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']
	vData = pd.read_table('supercomp.psl',sep='\t',header = None,names=Header_names,skiprows=5)

	#Make list of transcripts
	transcripts = np.unique(list(vData['Q name']))
	transcript_lengths = np.unique(list(vData['Q size']))

	#Metric 1 - If each transcript has a one to one mapping with ST then we should have as many lines as transcripts
	mapping = 0
	if(len(transcripts) == len(vData)): mapping =1

	#Metric 2 - The fraction fo each transcript mapped (i.e the number of bases)
	fraction = {} # Dictonary of fraction mapped
	for i in range(0,len(vData)):		
		if vData.iloc[i,9] in fraction:  fraction[vData.iloc[i,9]] += (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10]) #If key already in dictionary sum fractions
		fraction[vData.iloc[i,9]] = (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10])

	#Metric 3 - How compactified is the super transcript compared to reference
	tot_trans = 0
	for lengo in transcript_lengths:
		tot_trans += int(lengo)
	if(len(vData) > 0) : ST_len = int(vData.iloc[0,14])
	else: 
		ST_len = 1
		tot_trans = 1	
	compact = ST_len/tot_trans


	#Get Annotation of transcripts on SuperTranscript from BLAT
	anno = ""
	for j in range(0,len(vData)):
                tStarts = vData.iloc[j,20].split(",")
                blocksizes  = vData.iloc[j,18].split(",")
                for k in range(0,len(tStarts)-1): #Split qStarts, last in list is simply blank space
                        anno = anno + gene_name + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(int(tStarts[k]) +1) + '\t' + str(int(tStarts[k]) + int(blocksizes[k]) + 1) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + 'gene_id "' + gene_name +'"; transcript_id "' + vData.iloc[j,9] + '";'   + '\n'
		
	return(mapping, fraction, anno, compact)

if __name__=='__main__':
	if(len(sys.argv) != 2):
		print("Checker function requires one argument")
		print("The genome whose super transcripts you wish to quality check")
	else:
		Checker(sys.argv[1])
