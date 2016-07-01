#Visualise a given gene in your super transcript
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
import sys
import os
from matplotlib.pyplot import cm
from matplotlib import gridspec

##################################################
###### Visualise blocks in SuperTranscript #######
##################################################

def Visualise(gene_name):
	
	print("Producing Super Files\n")
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

	#Read in GFF file
	gff_data = pd.read_table('SuperDuper.gff',sep='\t',header=None)

	#Subset GFF just for gene
	gff_data = gff_data.loc[gff_data.ix[:,0] == gene_name,:]

	#Make list of transcripts
	transcripts = np.unique(list(vData['Q name']))

	#First pass attempt, make a stacked bar chart of all the nodes in C in order, with size dependent on length of string, iterate through coloirs

	#Get Super Transcript Length
	ST_length = int(vData.iloc[0,14])

	#SuperTranscript as one block
	#plt.barh(len(transcripts),ST_length,color='Black',left=0)

	gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
	ax1=plt.subplot(gs[0])
	accum = 0
	for row in range(0,len(gff_data)):
		size = int(gff_data.iloc[row,4]) - int(gff_data.iloc[row,3])
		plt.barh(len(transcripts),size,color='#ffc024',left=accum,alpha=0.8)
		accum=accum+size
		if(row > 0): plt.axvline(int(gff_data.iloc[row,3]),linestyle='dashed',color='black',linewidth=0.5)

	plot_dict = {}
	col_dict = {}
	labs = []

	col2=iter(cm.plasma(np.linspace(0,1,len(transcripts))))
	for i,key in enumerate(transcripts):
        	plot_dict[key] = i
        	col_dict[key] = next(col2)
        	lab = "T"+str(i)
        	labs.append(lab)

	#Empty vector to store coverage
	coverage = np.zeros(ST_length)

	for irow in range(0,len(vData)):

        	#Get blocks
		block_sizes = (vData.iloc[irow,18]).rstrip(',').split(',')
		tStarts = (vData.iloc[irow,20]).rstrip(',').split(',')
		qName = vData.iloc[irow,9]

        	#Print transcripts blocks
		for block in range(0,len(block_sizes)):
			si = int(block_sizes[block])
			left = int(tStarts[block])
			plt.barh(plot_dict[qName],si,color=col_dict[qName],left=left,alpha=0.7)
			#Sum up coverage
			for i in range(left,left+si):
				coverage[i] = coverage[i] + 1

	plt.setp(ax1.get_xticklabels(), visible=False)

	ind=np.arange(len(transcripts)+1)
	width=0.8
	labs.append('Super')
	plt.yticks(ind + width/2.,labs)
	plt.title('Breakdown of Super Transcript')
	plt.ylabel('Transcripts')


	#################################
	# Coverage Histogram Underneath #
	#################################

	#For a super block, how many transcripts cover that area....
	ax2=plt.subplot(gs[1],sharex=ax1)
	x = np.arange(ST_length)
	plt.bar(x,coverage,alpha=0.7)
	plt.xlim([0,ST_length+1])

	#Fix x-axes
	ax2.set_yticklabels([])
	plt.xlabel('Bases')
	plt.ylabel('Coverage')
	plt.savefig("Visualise.pdf")
	plt.show()

	###########################################
	# Create GFF with transcript annotation ###
	###########################################

	fg = open(gene_name + ".gff","w")
	#For each blat alignment
	for j in range(0,len(vData)): 
		qStarts = vData.iloc[j,19].split(",")
		blocksizes  = vData.iloc[j,18].split(",")
		for k in range(0,len(qStarts)-1): #Split qStarts, last in list is simply blank space
			fg.write(gene_name + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + qStarts[k] + '\t' + str(int(qStarts[k]) + int(blocksizes[k])) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + 'gene_id "' + gene_name +'"; transcript_id "' + vData.iloc[j,9] + '";'   + '\n')
	fg.close()
				


if __name__=='__main__':
	#Make argument parser
        parser = argparse.ArgumentParser()

        #Add Arguments
        parser.add_argument("GeneName",help="The name of the gene whom you wish to view")
        args= parser.parse_args()
        if(not os.path.isfile(args.GeneName + ".fasta")):
                print("No fasta file for gene/cluster of interest\n")
                sys.exit()
        if(not os.path.isfile("SuperDuper.fasta")):
                print("No fasta file for SuperTranscript\n")
                sys.exit()
        if(not os.path.isfile("SuperDuper.gff")):
                print("No annotation file for SuperTranscript\n")
                sys.exit()

        Visualise(args.GeneName)
