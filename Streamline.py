#A quick little first pass script to read a BLAT output and then:
# 1) Determine block sequences
# 2) Construct graph structure that stores each block with eges detailing how blocks are connect within transcripts
# 3) Sort blocks from the graph into topological order
# 4) Read sequence for each block to give the SuperTranscript!

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import time
import numpy as np
import sys

sys.setrecursionlimit(10000) # 10000 is an example, try with different values


#Start Clock for timing
start_time = time.time()

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


#################################
# Make test strings for debugging
#################################
#block_seq = [ ['T','A','C','G'], ['C','A','C','T'], ['G','G','C','T'], ['A','C'], ['C','T'], ['C']]
#tName = ['T1','T2','T3','T2','T3','T3']
#qName = ['T1','T2','T3','T1','T2','T1']
#tStart = [0,0,0,1,2,2]
#qStart = [0,0,0,1,2,2]
#transcripts = {'T1':'TACG','T2':'CACT','T3':'GGCT'}

#####################
# Construct the Graph 
######################
G= nx.DiGraph()
Node_index=0
node_dict = {} #A place to find the index of a node
for key in transcripts:
	node_dict[key] = np.full(len(list(transcripts[key])),-1,dtype=np.int64)
	
pn_id = 0 #Previous Node Id

#############################
# ADD NODES AND EDGES #######
#############################

#Add a node for every base in each transcript
for key in transcripts:
	for i in range(0,len(transcripts[key])):
		#Add node to graph
		G.add_node(Node_index)

              	#Add attributes
		G.node[Node_index]['Base'] = transcripts[key][i]

		#Add coordinates
		for k in transcripts:
			if(k == key): G.node[Node_index][k] = i
			else: G.node[Node_index][k] = None

		#Add to dictionary
		node_dict[key][i] = Node_index

		Node_index=Node_index + 1


#Now we want to merge the nodes in blocks
for i in range(0,len(block_seq)): #Loop through every block sequence

	if(tName[i] == qName[i]): continue #If comparing transcript to itself the nodes are already made

	for j in range(0,len(block_seq[i])): #Loop through every base in a given block 

		#print(nx.nodes(G))
	
		#labels=dict((n,d['Base']) for n,d in G.nodes(data=True))
		#nx.draw(G,labels=labels)
		#plt.show()

		#Co-ordinate for base in transcript and query sequences
		tpos = j + tStart[i]
		qpos = j + qStart[i]
	
		#Node index for q base and t base
		qnid = node_dict[qName[i]][qpos] #Node index in graph for base in query
		tnid = node_dict[tName[i]][tpos] #Node index in graph for base in transcript
		

		#Check if node already there either with the transcript id or query id

		#print("### MERGE START ###")
		#print("Started with:",qnid)
		#print("Merging Node: ", G.node[qnid])
		#print("To: ", G.node[tnid])
		#print("Which is:",tnid)	

		#If they are not the same node, we need to merge them and add the same edges, redirect the query node to the transcript node
		if(qnid != tnid): 

			#Redirect incoming edges
			for n1,n2 in G.in_edges([qnid]): #For each pair of nodes where there is an edge for query nodes 
				G.add_edge(n1,tnid)

			#Redirect Outgoing edges
			for n1,n2 in G.out_edges([qnid]): #For each pair of nodes where there is an edge for query nodes
				G.add_edge(tnid,n2)				

			#Merge attributes, without overwriting transcript node, first for query node
			G.node[tnid][qName[i]] = qpos	
			
			#Loop through attributes in query node and add if not none
			for key in transcripts:
				if(key == qName[i]):
					continue

				#Only override transcript attributes if not empty
				if(G.node[qnid][key] is not None):
					if(key != tName[i]): G.node[tnid][key] = G.node[qnid][key] #Don't replace transcript position
					#Update Dictionary
					node_dict[key][G.node[qnid][key]] = tnid
				
					
				
			#################
			# Remove old node
			#################

			#Remove query node since we have no merged it to the transcript node
			G.remove_node(qnid)

			#Change Dictionary Call for query node
			node_dict[qName[i]][qpos] = tnid



###################################################
## Add Edges between adjacent transcript nodes ####
###################################################
for key in node_dict:
        for j in range(0,len(node_dict[key])-1):
                G.add_edge(node_dict[key][j],node_dict[key][j+1])





#Manually break a cycle
#Check for cycles
#print(len(list(nx.simple_cycles(G))))
cycles = list(nx.simple_cycles(G)))
#G.remove_node(666)

########
#Check we can get transcripts back from graph
########

for key in transcripts:
	path_alternatives = []
	#print("Transcript: ", transcripts[key])
	pathl = nx.all_simple_paths(G,node_dict[key][0],node_dict[key][-1])
	for path in pathl:
		seq= ''
		#print(path)
		for b in path:
			#print(G.node[b]['Base'])
			seq = seq + G.node[b]['Base']
		#print(seq)
		path_alternatives.append(seq)

	#Loop through and check the strings match
	for path in path_alternatives:
		if(path == transcripts[key]): print("Found path for: ", key)
		



#Order base in graph
#base_order = nx.topological_sort_recursive(G)
#seq =''
#for index in base_order:
#	seq = seq + G.node[index]['Base']
#print(seq)	

#For Debugging
#labels=dict((n,d['Base']) for n,d in G.nodes(data=True))
#nx.draw(G,labels=labels)
#nx.draw_circular(G,labels=labels)
#nx.draw_spectral(G,labels=labels)
#plt.show()
	
#print(len(G))
#print("Finding Node")
#find_node = [attrdict for n,attrdict in G.node.items() if attrdict[qName[1]] == qStart[1] ]
#print(find_node)
#print(qName[1])
#print(len(find_node))

#for i in range(0,12):
#	print(G.node[i])				

print("---- %s seconds ----" %(time.time()-start_time))
