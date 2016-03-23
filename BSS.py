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

#################################
# Make test strings for debugging
#################################
block_seq = [ ['T','A','C','G'], ['C','A','C','T'], ['G','G','C','T'], ['A','C'], ['C','T'], ['C']]
tName = ['T1','T2','T3','T2','T3','T3']
qName = ['T1','T2','T3','T1','T2','T1']
tStart = [0,0,0,1,2,2]
qStart = [0,0,0,1,2,2]
transcripts = ['T1','T2','T3']

#####################
# Construct the Graph 
######################
G= nx.DiGraph()
Node_index=0
node_dict = {} #A place to find the index of a node
pn_id = 0 #Previous Node Id

#############################
# ADD NODES #################
#############################

#Add nodes, each node is a base
for i in range(0,len(block_seq)): #Loop through every block sequence
	for j in range(0,len(block_seq[i])): #Loop through every base in a given block (N.B. One of the blocks is the whole transcript)

		#print(nx.nodes(G))
	
		#labels=dict((n,d['Base']) for n,d in G.nodes(data=True))
		#nx.draw(G,labels=labels)
		#plt.show()

		dict_name_t = tName[i] + ':' + str(j + tStart[i])
		dict_name_q = qName[i] + ':' + str(j + qStart[i])		
	

		#Check if node already there either with the transcript id or query id
		if(node_dict.get(dict_name_t) and node_dict.get(dict_name_q)): #If already have a node for both bases

			#print("### MERGE START ###")
			#print("Started with:",node_dict.get(dict_name_q))
			#print("Transcript Node: ", dict_name_t)
			#print("Query Node: ", dict_name_q)
			#print("Merging Node: ", G.node[node_dict[dict_name_q]])
			#print("To: ", G.node[node_dict[dict_name_t]])	

			#If they are not the same node, we need to merge them and add the same edges, redirect the query node to the transcript node
			if(node_dict[dict_name_t] != node_dict[dict_name_q]): 

				#For each pair of nodes where there is an edge
				#print(G.in_edges([node_dict[dict_name_q]]))
				#print(G.out_edges([node_dict[dict_name_q]]))

				#Redirect incoming edges
				for n1,n2 in G.in_edges([node_dict[dict_name_q]]): #For each pair of nodes where there is an edge for query nodes 
					G.add_edge(n1,node_dict[dict_name_t])
				#Redirect Outgoing edges
				for n1,n2 in G.out_edges([node_dict[dict_name_q]]): #For each pair of nodes where there is an edge for query nodes
					G.add_edge(node_dict[dict_name_t],n2)				

				#Merge attributes, without overwriting transcript node
				G.node[node_dict[dict_name_t]][qName[i]] = j+ qStart[i]	
			
				#Loop through attributes in query node and add if not none
				for key in transcripts:
					if(key == qName[i] or key == tName[i]):
						continue

					#Only override transcript attributes if not empty
					if(G.node[node_dict[dict_name_q]][key] is not None):
						G.node[node_dict[dict_name_t]][key] = G.node[node_dict[dict_name_q]][key]

						#Update Dictionary
						dict_name_temp = key + ":" + str(G.node[node_dict[dict_name_q]][key])
						node_dict[dict_name_temp] = node_dict[dict_name_t]
	
				
				#print ("Removing node:",  node_dict[dict_name_q]) 
				#print ("Merged to: ", G.node[node_dict[dict_name_t]])
				#print ("Which is node: ", node_dict[dict_name_t])
				
				#Remove query node since we have no merged it to the transcript node
				G.remove_node(node_dict[dict_name_q])

				#Change Dictionary Call for query node
				node_dict[dict_name_q] = node_dict[dict_name_t]


			#print("Changed to:", node_dict.get(dict_name_q))
			
			
			pn_id = node_dict[dict_name_t]
			continue

		elif(node_dict.get(dict_name_t)): #If there is already a node for that transcript position, append query info
			print("Already transcript node")
			temp_id = node_dict[dict_name_t]
			G.node[temp_id][qName[i]] = j+ qStart[i]
			pn_id = temp_id
			node_dict[dict_name_q] = temp_id

		elif(node_dict.get(dict_name_q)): #Already a node for that query transcript position, append query info
			print("Already query node")
			temp_id = node_dict[dict_name_q]
			G.node[temp_id][tName[i]] = j+ tStart[i]
			pn_id = node_dict[dict_name_q]
			node_dict[dict_name_t] = temp_id

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
			node_dict[dict_name_t] = Node_index
			node_dict[dict_name_q] = Node_index

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
nx.draw(G,labels=labels)
plt.show()
	
#print(len(G))
#print("Finding Node")
#find_node = [attrdict for n,attrdict in G.node.items() if attrdict[qName[1]] == qStart[1] ]
#print(find_node)
#print(qName[1])
#print(len(find_node))

#for i in range(0,12):
#	print(G.node[i])				
