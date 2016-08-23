#A script to systematically check for each gene whether the SuperTranscript builder worked ok
import multiprocessing, logging
from multiprocessing import Pool
from multiprocessing import Process
import argparse
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.pyplot import cm
plt.style.use('seaborn-deep')
from matplotlib import gridspec
import pickle
import time
from matplotlib.backends.backend_pdf import PdfPages

################################################################
###### Visualise blocks and metricise in SuperTranscript #######
################################################################

def Checker(genome,ncore):
	start_time = time.time()
	print("Finding list of genes")
	genes=[]
	whirls=[]
	f = open(genome,'r')
	for line in f:
		if('>' in line):
			genes.append((line.split('>')[1]).split(" ")[0])
			whirls.append(int((line.split('>')[1]).split(":")[-1]))

	# BY POOL
	pool = Pool(processes=ncore)
	result = pool.map_async(FindMetrics,genes)
	pool.close()
	pool.join()
	results = result.get()

	metrics = {}
	for i,gene in enumerate(genes):
		mapping,fraction,anno,compact,ntran = results[i]
		metrics[gene] = [mapping,fraction,anno,compact,ntran]

	#Fraction of genes where we get a one to one mapping
	mapp_frac = []
	frac_covered = []
	compactify = []
	trans = []
	for key in metrics:
		mapp_frac.append(metrics[key][0])
		compactify.append(metrics[key][3])
		trans.append(metrics[key][4])	
		for key2 in metrics[key][1]:
			frac_covered.append(metrics[key][1][key2])

	#mapp_frac = mapp_frac / len(genes)
	frac_covered = np.asarray(frac_covered)

	#Calculate the number of clusters with a one to one mapping with >1 transcript
	multiclust = 0
	mapp_mult=0
	for i in range(0,len(trans)):
		if(trans[i] > 1): #More than 1 contig
			multiclust += 1
			mapp_mult += mapp_frac[i]
		
	if(multiclust>0): mapp_mult=mapp_mult/multiclust

	with PdfPages('LogOut.pdf') as pdf:

		#Now lets pring some info
		p1 = plt.barh(0,mapp_mult,alpha=0.75,label='Fully Mapped',color='#55A868')
		p2 = plt.barh(0,(1.-mapp_mult),alpha=0.75,left=mapp_mult,label='Partially Mapped',color='#4C72B0')
		plt.yticks([],[])
		plt.title('Contiguous Mapping of transcripts to SuperTranscript per Cluster')
		plt.legend(loc='best',frameon=False)
		plt.xlim(0.,1.)
		plt.ylim(0,0.5)
		pdf.savefig()
		plt.close()

		#Visualise the number of transcripts per gene
		tt = np.asarray(trans)
		cbins= [1,2,3,4,5,6,7,8,9,10,20,500]
		h = np.histogram(tt,cbins)
		plt.bar(range(len(cbins)-1),h[0],alpha=0.75,width=1,color='#C44E52')
		xlab=['1','2','3','4','5','6','7','8','9','10','11-20','21-500']
		plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5],xlab)
		plt.xlabel("Number of transcripts in cluster")
		plt.ylabel("Frequency")
		plt.title('Distribution of transcripts per cluster')
		pdf.savefig()
		plt.close()

		#Plot the distribution for metric 2
		n, bins, patches = plt.hist(frac_covered,20,alpha=0.75,range=[0,1],color='#8172B2')
		plt.xlabel("Fraction of transcript covered in SuperTranscript")
		plt.ylabel("Frequency")
		plt.title('Fraction of transcripts encroporated into SuperTranscript')
		pdf.savefig()
		plt.close()

		#Plot Distribution for Metric 3
		compactified = np.asarray(compactify)
		#print(compactified)
		n,bins,patches = plt.hist(compactified,20,alpha=0.75,range=[0,1],color='#4C72B0')
		plt.xlabel("Super Transcript Length/Sum of bases in Transcripts")
		plt.ylabel("Frequency")
		plt.title("Compactifcation of transcripts in SuperTranscript")
		pdf.savefig()
		plt.close()

		#Save also the whirl info
		whirls = np.asarray(whirls)
		n,bins,patches = plt.hist(whirls,10,alpha=0.75,range=[0,10],color='#8172B2')
		plt.xlabel("Number of whirls in cluster")
		plt.ylabel("Frequency")
		plt.title("General Whirl Statistics")
		pdf.savefig()
		plt.close()


		# We can also set the file's metadata via the PdfPages object:
		d = pdf.infodict()
		d['Title'] = 'Metrics for the SuperTranscript Build'
		d['Author'] = u'Anthony Hawkins'
		d['ModDate'] = datetime.datetime.today()

	#Save the metrics as pickle file
	pickle.dump(metrics,open("Metrics.pkl","wb"))

	#Now lets save the annotation to file
	fg = open("SuperDuperTrans.gff","w")
	for key in metrics:
		fg.write(metrics[key][2])
	fg.close()
	print("ANNOCHECKED ---- %s seconds ----" %(time.time()-start_time))
	
def FindMetrics(gene_name):

	#EXTRACT GENE FROM SUPER
	gene_string=""
	#Find gene in genome
	f= open("SuperDuper.fasta","r")
	for line in f:
		if('>' in line):
			if( gene_name == (line.split('>')[1]).split(" ")[0]):
				gene_string=next(f)
				break
	f.close()

	fgf = "Super_" + str(gene_name) + ".fasta"
	fs= open(fgf,"w") 
	fs.write(">" + gene_name + "\n")
	fs.write(gene_string)
	fs.close()


	#Match transcripts to super transcript
	print("Producing match to super transcript")
	BLAT_command = "blat Super_%s.fasta %s.fasta super_%s.psl" %(gene_name,gene_name,gene_name)
	os.system(BLAT_command)

	#First read in BLAT output:
	Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']
	fgp = 'super_'+str(gene_name)+'.psl'
	vData = pd.read_table(fgp,sep='\t',header = None,names=Header_names,skiprows=5) 

	#Make list of transcripts
	transcripts = np.unique(list(vData['Q name']))
	transcript_lengths = np.unique(list(vData['Q size']))
	ntran = len(transcripts)

	#Metric 1 - Check for each transcript whether or not there is a row in the psl where the transcript maps at least 95% to the SuperTranscript
	mapping = 0
	tran_map = 0

	#Metric 2 - The fraction fo each transcript mapped (i.e the number of bases)
	fraction = {} # Dictonary of fraction mapped
	for i in range(0,len(vData)):		
		if vData.iloc[i,9] in fraction:  fraction[vData.iloc[i,9]] += (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10]) #If key already in dictionary sum fractions
		fraction[vData.iloc[i,9]] = (int(vData.iloc[i,12]) - int(vData.iloc[i,11]))/int(vData.iloc[i,10])
		if(vData.iloc[i,0] >= 0.95*vData.iloc[i,10]): tran_map=tran_map+1 #Metric1

	if(tran_map == ntran): mapping =1

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
                        anno = anno + gene_name + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(int(tStarts[k]) +1) + '\t' + str(int(tStarts[k]) + int(blocksizes[k])) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + 'gene_id "' + gene_name +'"; transcript_id "' + vData.iloc[j,9] + '";'   + '\n'

	#Remove the psl and fasta file needed for the blat
	dcom='rm Super_%s.fasta' %(gene_name)
	os.system(dcom)
	dcom = 'rm super_%s.psl' %(gene_name)
	os.system(dcom)
		
	return(mapping, fraction, anno, compact,ntran)

if __name__=='__main__':
	#Make argument parser
	parser = argparse.ArgumentParser()

        #Add Arguments
	parser.add_argument("SuperFile",help="The name of the SuperDuper.fasta file created by SuperTranscript")
	parser.add_argument("--cores",help="The number of cores you wish to run the job on (default = 1)",default=1,type=int)
	args= parser.parse_args()
	Checker(args.SuperFile,args.cores)
