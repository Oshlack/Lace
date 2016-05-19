#A cute script to visualise the super transcripts blocks and the blocks from the transcripts which build it
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

#Match transcripts to super transcript
print("Producing match to super transcript")
BLAT_command = "./blat Super.fasta CLC_stranded.fasta supercomp.psl"
os.system(BLAT_command)

#First read in BLAT output:
Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']
vData = pd.read_table('supercomp.psl',sep='\t',header = None,names=Header_names,skiprows=5)

#Read in GFF file
gff_data = pd.read_table('Super.gff',sep='\t',header=None)

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
plt.bar(x,coverage,color='slategray',alpha=0.7)
plt.xlim([0,ST_length+1])

#Fix x-axes
ax2.set_yticklabels([])
plt.xlabel('Bases')
plt.ylabel('Coverage')
plt.show()



