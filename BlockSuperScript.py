#A quick little first pass script to read a BLAT output and then:
# 1) Determine block sequences
# 2) Construct graph structure that stores each block with eges detailing how blocks are connect within transcripts
# 3) Sort blocks from the graph into topological order
# 4) Read sequence for each block to give the SuperTranscript!

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#####################
#Read in transcripts
#####################

fT = open('ESR1_transcripts_for_Michael.fasta','r')
transcripts = {}
tName = ''

for line in fT:
	if(">" in line): #Name of 
		tName = line.split('\n')[0].split('\r')[0].lstrip('>')
		transcripts[tName] = ''
	
	else:
		transcripts[tName] = transcripts[tName] + line.split('\n')[0].split('\r')[0]


###################################
# Loop pairwise through transcripts
###################################

#This may well be changed/ skipped/ parallelised
#BLAT_command = "blat ESR1_transcripts_for_Michael.fasta ESR1_transcripts_for_Michael.fasta outall.psl"
#os.system(BLAT_command)


#First read in BLAT output:
Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']

bData = pd.read_table('outall.psl',sep='\t',header = None,names=Header_names,skiprows=5)
print(bData)

###############
#Now extract the sequences from the blocks using the coordinates from BLAT
###############

f = open('ENST00000338799.5','r')
sequence = ''
for i, line in enumerate(f):
	if(i>0): sequence = sequence + line.split('\n')[0].split('\r')[0]

seq = list(sequence) # Transform string into list of characters (bases)

block_seq= []
tName = [] #The name of the transcript to which the block is shared in 
tStart = [] #Start co-ordinate of the block in the transcript
qStart = []

for i in range(0,len(bData)):

	block_sizes = (bData.iloc[i,18]).rstrip(',').split(',')
	qStarts = (bData.iloc[i,19]).rstrip(',').split(',')
	tStarts = (bData.iloc[i,20]).rstrip(',').split(',')

	for j in range(0,len(qStarts)):
		block_seq.append(seq[int(qStarts[j]):(int(qStarts[j])+int(block_sizes[j]))])
		tName.append(bData.iloc[i,13])
		tStart.append(int(tStarts[j]))
		qStart.append(int(qStarts[j]))

#print(bData)

#####################
# Construct the Graph 
######################

#Add nodes, each node is a block
G = nx.DiGraph()
for i in range(0,len(block_seq)): #Loop through every block sequence

		#Add node to graph
		G.add_node(i)

		#Add attributes
		G.node[i]['Transcript'] = tName[i] #Which transcript is in the block
		

				
#Add edges - for each block pair, if the start of one block is after the other then we say that they sequentially follow the other
for i in range(0,len(block_seq)):
	for j in range(0,len(block_seq)):
		if(i == j): continue
		if(qStart[j] > qStart[i]):	G.add_edge(i,j)



#nx.draw(G)
#plt.show()	

#Order blocks in graph
block_order = nx.topological_sort(G)
print(block_order)		
		
				
