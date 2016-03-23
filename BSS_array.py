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

#Add nodes, each node is a base
for i in range(0,len(block_seq)): #Loop through every block sequence
	for j in range(0,len(block_seq[i])): #Loop through every base in a given block (N.B. One of the blocks is the whole transcript)

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
		
		#print("Node index of mystery: ", node_dict['ENST00000443427.1'][1356])

		#Check if node already there either with the transcript id or query id
		if(tnid > -1 and qnid > -1): #If already have a node for both bases

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
			
				if(qnid == 977):
					print ("Removing node:",  qnid)
					print(G.node[qnid])
					print("Replaced by:", tnid)
					print(G.node[tnid])
					
					
				#Loop through attributes in query node and add if not none
				for key in transcripts:
					if(key == tName[i] or key == qName[i]):
						continue

					#Only override transcript attributes if not empty
					#print(nx.nodes(G))
					print("Node id ", qnid)
					print("Key ", key)
					print("Q pos", qpos)
					print("Name",qName[i])
					#print(G.node[qnid][key])

					if(G.node[qnid][key] is not None):

						print("TName:",tName[i])
						print("QName:",qName[i])
						print("Key: ", key)
						print("Pos: ", G.node[qnid][key])
						print("Node Index ", node_dict[key][G.node[qnid][key]])

						G.node[tnid][key] = G.node[qnid][key]

						#Update Dictionary
					#	print(G.node[qnid][key])
						node_dict[key][G.node[qnid][key]] = tnid
					#	print(node_dict[key][G.node[qnid][key]])
				
					
				#print ("Merged to: ", G.node[tnid])
				
				#################
				# Remove old node
				#################

				#Remove query node since we have no merged it to the transcript node
				G.remove_node(qnid)

				#Change Dictionary Call for query node
				node_dict[qName[i]][qpos] = tnid

				if(qnid == 977):
					print(tnid)
					print ("Should Now read:", node_dict[qName[i]][qpos])
					#Check all dictionary postions
					for key in transcripts:
						print("Key: ", key)
						print("Pos: ", G.node[tnid][key])
						if(G.node[tnid][key] is not None): print("Node Id:", node_dict[key][G.node[tnid][key]])


			#print("Changed to:", node_dict[qName[i]][qpos])
			
			
			pn_id = tnid
			continue

		elif(tnid > -1): #If there is already a node for that transcript position but not for query
			G.node[tnid][qName[i]] = qpos #Update attribute on existing node
			pn_id = tnid #Update previous node
			node_dict[qName[i]][qpos] = tnid #Update previous dictionary

		elif(qnid > -1): #If there is already a node for that transcript position but not for query
                        G.node[qnid][tName[i]] = tpos
                        pn_id = qnid
                        node_dict[tName[i]][tpos] = qnid

		else: #Add a new node

			#Add node to graph		
			G.add_node(Node_index)

			#Add attributes
			G.node[Node_index]['Base'] = block_seq[i][j]

			for key in transcripts:
				if(key == tName[i]): G.node[Node_index][key] = j + tStart[i] 		
				elif(key == qName[i]): G.node[Node_index][key] = j + qStart[i]
				else: G.node[Node_index][key] = None


			
			#Add Edge now to previous node in query sequence
			if(j > 0):
				G.add_edge(pn_id,Node_index)
	
			pn_id = Node_index
			node_dict[tName[i]][tpos] = Node_index
			node_dict[qName[i]][qpos] = Node_index

			#Add one to node iterator, if we have actually added a node
			Node_index = Node_index + 1




#Order base in graph
base_order = nx.topological_sort_recursive(G)
seq =''
for index in base_order:
	seq = seq + G.node[index]['Base']
print(seq)	


#For Debugging
labels=dict((n,d['Base']) for n,d in G.nodes(data=True))
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
