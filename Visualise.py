#A cute script to visualise the super transcripts blocks and the blocks from the transcripts which build it
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.pyplot import cm
import seaborn as sns
#plt.style.use('fivethirtyeight')

##################################################
###### Visualise blocks in SuperTranscript #######
##################################################

#Load graph
C = nx.read_gpickle('Simples.pkl')

#Match transcripts to super transcript
if(os.path.isfile('supercomp.psl') == False):
	BLAT_command = "./blat Super.fasta ESR1_transcripts_for_Michael.fasta -maxGap=0 -minIdentity=100 -maxIntron=0  supercomp.psl"
	os.system(BLAT_command)

#First read in BLAT output:
Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']
vData = pd.read_table('supercomp.psl',sep='\t',header = None,names=Header_names,skiprows=5)

#Make list of transcripts
transcripts = np.unique(list(vData['Q name']))

#First pass attempt, make a stacked bar chart of all the nodes in C in order, with size dependent on length of string, iterate through coloirs
N = nx.number_of_nodes(C)
node_sizes = []
colour=iter(cm.viridis(np.linspace(0,1,N)))

left =0
for n in C.nodes():
        base_string = list(C.node[n]['Base'])
        node_sizes.append(len(base_string))
        c = next(colour)
        plt.barh(len(transcripts),len(base_string),color=c,left=left)
        left = left + len(base_string)

plot_dict = {}
col_dict = {}
labs = []

col2=iter(cm.plasma(np.linspace(0,1,len(transcripts))))
for i,key in enumerate(transcripts):
        plot_dict[key] = i
        col_dict[key] = next(col2)
        lab = "T"+str(i)
        labs.append(lab)

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

#Fix y-axis
ind=np.arange(len(transcripts)+1)
width=0.8
labs.append('Super')
plt.yticks(ind + width/2.,labs)
plt.title('Breakdown of Super Transcript')
plt.xlabel('Bases')
plt.ylabel('Transcripts')
#plt.figure(facecolor='white')
plt.show()
