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


###################################################
# Loop pairwise through transcripts and BLAT allign
###################################################

#This may well be changed/ skipped/ parallelised
#BLAT_command = "blat ESR1_transcripts_for_Michael.fasta ESR1_transcripts_for_Michael.fasta outall.psl"
#os.system(BLAT_command)


#First read in BLAT output:
Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']

bData = pd.read_table('outall.psl',sep='\t',header = None,names=Header_names,skiprows=5)
#print(bData)

###############
#Now extract the sequences from the blocks using the coordinates from BLAT
###############
block_seq= []
tName = [] #The name of the transcript to which the block is shared in 
qName = []
tStart = [] #Start co-ordinate of the block in the transcript
qStart = [] #Start co-ordinate of the block in query

for i in range(0,len(bData)):

	seq=list(transcripts[bData.iloc[i,9]]) #Get sequence from query name
	block_sizes = (bData.iloc[i,18]).rstrip(',').split(',')
	qStarts = (bData.iloc[i,19]).rstrip(',').split(',')
	tStarts = (bData.iloc[i,20]).rstrip(',').split(',')

	for j in range(0,len(qStarts)):
		block_seq.append(seq[int(qStarts[j]):(int(qStarts[j])+int(block_sizes[j]))])
		tName.append(bData.iloc[i,13])
		tStart.append(int(tStarts[j]))
		qStart.append(int(qStarts[j]))
		qName.append(bData.iloc[i,9])

#####################
# Construct the Graph 
######################
G= nx.DiGraph()
Node_index=0
to_add_edges = []
node_dict = {} #A place to find the index of a node

#############################
# ADD NODES #################
#############################

#Add nodes, each node is a base
for i in range(0,len(block_seq)): #Loop through every block sequence
	for j in range(0,len(block_seq[i])): #Loop through every base in a given block (N.B. One of the blocks is the whole transcript)
	
		dict_name_t = tName[i] + ':' + str(j + tStart[i])
		dict_name_q = qName[i] + ':' + str(j + qStart[i])		



		#Add node to graph		
		G.add_node(Node_index)

		#Add attributes
		G.node[Node_index]['Base'] = block_seq[i][j]

		for key in transcripts:
			if(key == tName[i]): G.node[Node_index][key] = j + tStart[i] 		
			elif(key == qName[i]): G.node[Node_index][key] = j + qStart[i]
			else: G.node[Node_index][key] = None

		#Add Edge now to previous node in query sequence
		#if(Node_index > 0): G.add_edge(Node_index-1,Node_index)
	
		#if(qName[i] != tName[i]): #If query maps to a transcript other than itself then need to add edge
		#	to_add_edges.append(Node_index)
				

		#Add one to node iterator
		Node_index = Node_index + 1

print(len(to_add_edges))		

#Loop through the nodes with edges to add		
#for node in to_add_edges:
#	for key in Transcripts:
#		if(G[node][key] != None) : #For the nodes
			


#nx.draw(G)
#plt.show()	

#Order blocks in graph
#block_order = nx.topological_sort(G)
#print(block_order)		
		
#print(len(G))
print("Finding Node")
find_node = [attrdict for n,attrdict in G.node.items() if attrdict[qName[1]] == qStart[1] ]
print(find_node)
print(qName[1])
print(len(find_node))

#for i in range(0,6):
#	print(G.node[i])				
