#A little script to construct an annotation from an SJ.out.tab file

import pandas as pd
import os, sys


#First read in the SJ.out.tab as sys.argv[1]
sj = pd.read_csv(sys.argv[1],sep='\t',header=None,names=['Gene','Start','End','strand','intron motif','Annotated','Unique','Multi','Overhang'])
print(sj.head())

#Read in SuperTranscript fasta file
sf = open(sys.argv[2],'r')
glength = {} #A dictionary holding the SuperTranscript gene lengths
gene=''
for line in sf:
	if('>' in line):
		gene= (line.split(' ')[0]).split('>')[1]
		glength[gene] = ''
	else:
		glength[gene] = glength[gene] + line.split('\n')[0].split('\r')[0]

#Create gtf file
gtf = open('Spliced.gtf','w')
curr_gene = sj.iloc[0,0]	
slist = {}


#Make a dictionary for each gene, which holds a list of splice junction start and end points
#For each row
for i in range(0,len(sj['Gene'])):
	curr_gene = sj.iloc[i,0]

	if(curr_gene not in slist.keys()): slist[curr_gene] = [1]

	if((sj.iloc[i,7] + sj.iloc[i,8]) > 10): #More than 10 reads (either unique or multi spanninf junction)
		slist[curr_gene].append(int(sj.iloc[i,1]))
		slist[curr_gene].append(int(sj.iloc[i,2])+1) #This is the end of the exon-exon junction, so need to add one for the actual exon start place

#Now sort each list for each gene, only keep unique elements
for key in glength.keys():
	if(key in slist.keys()):
		slist[key] = list(set(slist[key]))
		slist[key].sort()

		#Now for each coord in list make a block
		for i in range(1,len(slist[key])):
			ann = str(key) + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(slist[key][i-1])  + '\t' + str(slist[key][i]-1) + '\t' + '.' + '\t' + '.' + '\t' + '0' + '\t' + '.' + '\n' #Note: need to minus one off end to account for the fact that the exon ends before the exon-exon boundary exists
			gtf.write(ann)


	#For the list splice junnction, we need to make a block from the last sj to the end of the ST
	if(key not in slist): last = 1
	else: last = slist[key][-1]
 
	ann = str(key) + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(last)  + '\t' + str(len(glength[key])) + '\t' + '.' + '\t' + '.' + '\t' + '0' + '\t' + '.' + '\n'

	gtf.write(ann)
