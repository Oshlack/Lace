#A quick little first pass script to read a BLAT output and then:
# 1) Determine block sequences
# 2) Construct graph structure that stores each block with eges detailing how blocks are connect within transcripts
# 3) Sort blocks from the graph into topological order
# 4) Read sequence for each block to give the SuperTranscript!

import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import time
import numpy as np
import sys
from matplotlib.pyplot import cm 
import os
sys.setrecursionlimit(10000) # 10000 is an example, try with different values


#A little function to loop recursively through linked nodes in Graph in order to replace the node stored due to merging/removal
#def Recursive_Replace(G,node_dict,transcripts,T,Q,tnid,id):
#	for key in transcripts:
#		if(G.node[id][key] is not None):
#                                        #Update Dictionary
#                                        node_dict[key][G.node[id][key]] = tnid
#					Recursive_Replace(G,node_dict,transcripts,T,Q,tnid,node_dict[key][G.node[id][key]])	
	



#Start Clock for timing
start_time = time.time()

#####################
#Read in transcripts
#####################

fT = open('ESR1_transcripts_for_Michael.fasta','r')
#fT = open('s.cerevisiae/Cluster-3231.0.fasta','r')
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
if(os.path.isfile('outall.psl') == False):
	BLAT_command = "./blat ESR1_transcripts_for_Michael.fasta ESR1_transcripts_for_Michael.fasta -maxGap=0 -minIdentity=100 -maxIntron=0 outall.psl"
	#BLAT_command = "./blat s.cerevisiae/Cluster-3231.0.fasta s.cerevisiae/Cluster-3231.0.fasta  -maxGap=0 -minIdentity=100 -maxIntron=0 outall.psl"
	os.system(BLAT_command)


#First read in BLAT output:
Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']

bData = pd.read_table('outall.psl',sep='\t',header = None,names=Header_names,skiprows=5)
#bData = pd.read_table('ESR1_exakt.psl',sep='\t',header = None,names=Header_names,skiprows=5)

###############
#Now extract the sequences from the blocks using the coordinates from BLAT
###############

block_seq= []
tName = [] #The name of the transcript to which the block is shared in 
qName = []
tStart = [] #Start co-ordinate of the block in the transcript
qStart = [] #Start co-ordinate of the block in query

pair_list = []

for i in range(0,len(bData)):
	#Check explicitly that there are no gaps
	if( bData.iloc[i,4] > 0 or bData.iloc[i,6] > 0):
		continue
	
	#Don't allign the transcripts against each other twice...
	#I.e. BLAT does T1 vs T2 but also T2 vs T1 (which should be the same give or take)
	
	#TName + QName
	paired = bData.iloc[i,13] + bData.iloc[i,9]
	if(paired in pair_list): continue
	else: pair_list.append(bData.iloc[i,9]+bData.iloc[i,13]) 

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

#Motherload
#block_seq = [['A','C','T'],['G','A'],['A','C','T'],['G','A'],['A','C','T'],['A','C','T'],['C','G','G'],['C','G','G'],['A','C','T'],['G','A'],['C','G','G'],['A','C','T'],['C','G','G'],['C','G','G']]
#tName = ['T1','T1','T1','T1','T1','T2','T2','T2','T2','T2','T2','T3','T3','T3']
#qName = ['T2','T2','T3','T3','T4','T4','T4','T4','T3','T3','T3','T4','T4','T4']
#tStart = [0,4,0,4,0,0,4,4,0,8,4,0,7,7]
#qStart = [0,8,0,4,0,0,4,8,0,4,7,0,4,8]
#transcripts = {'T1':'ACTGGAC','T2':'ACTGCGGAGA','T3':'ACTCGAGCGG','T4':'ACTTCGGACGG'}

#Example with a loop
#block_seq=[['A','C','T'],['A','C','T'],['A','C','T'],['A','C','T']]
#tName = ['T2','T2','T1','T1']
#qName = ['T1','T1','T2','T2']
#tStart = [0,0,0,4]
#qStart = [0,4,0,0]
#transcripts = {'T1':'ACTGACTA','T2':'ACTGCGCA'}
 


########################
# Construct the Graph ##
########################
G= nx.DiGraph()
#G=nx.MultiDiGraph()

Node_index=0
node_dict = {} #A place to find the index of a node
for key in transcripts:
	node_dict[key] = [-1] * len(list(transcripts[key]))

#############################
# ADD NODES FOR EACH BASE #######
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


###################################################
## Add Edges between adjacent transcript nodes ####
###################################################
for key in node_dict:
	for j in range(0,len(node_dict[key])-1):
		#Add a weight for the number of transcripts which make this edge, and also flag which transcripts are using this link
		wt= 1
		#If the edge already exists
		if (G.has_edge(node_dict[key][j],node_dict[key][j+1])): wt = G[node_dict[key][j]][node_dict[key][j+1]]['weight'] + 1
		G.add_edge(node_dict[key][j],node_dict[key][j+1],weight=wt,key=1)


#######################################
# Glue           ######################
#######################################

glue_nodes = 0
#Now we want to merge the nodes in blocks
for i in range(0,len(block_seq)): #Loop through every block sequence

	if(tName[i] == qName[i]): continue #If comparing transcript to itself the nodes are already made

	for j in range(0,len(block_seq[i])): #Loop through every base in a given block 

		#Co-ordinate for base in transcript and query sequences
		tpos = j + tStart[i]
		qpos = j + qStart[i]
	

		#Node index for q base and t base
		qnid = node_dict[qName[i]][qpos] #Node index in graph for base in query
		tnid = node_dict[tName[i]][tpos] #Node index in graph for base in transcript
		

		#Check if node already there either with the transcript id or query id


		if(qnid != tnid): #query and transcript node not the same

			#Consideration 1 - Whirls from repeated sections
			#Check if transcript node id already used for another base on the query string
			try:
				ll = node_dict[qName[i]].index(tnid)

			except ValueError:
    				ll = -1

			#Check whether the transcript node you are merging to, isnt already in the query string
			if(ll >= 0 and ll != qpos): continue		

			#If the node you are intending to merge is already merged to somewhere else on the transcript string, dont merge as can cause wirls
			if(qnid in node_dict[tName[i]]): continue 

			
			#print("### MERGE START ###")
			#print("Started with:",qnid)
			#print("Merging Node: ", G.node[qnid])
			#print("To: ", G.node[tnid])
			#print("Which is:",tnid)

			#Redirect incoming edges
			for n1,n2,d in G.in_edges([qnid],data=True): #For each pair of nodes where there is an edge for query nodes 
				G.add_edge(n1,tnid,d)

			#Redirect Outgoing edges
			for n1,n2,d in G.out_edges([qnid],data=True): #For each pair of nodes where there is an edge for query nodes
				G.add_edge(tnid,n2,d)				

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

			#Recursive check that no element in node dict contains the old node which is removed
			for key in node_dict:
				if(qnid in node_dict[key]):
					node_dict[key][node_dict[key].index(qnid)] = tnid


			glue_nodes = glue_nodes +1
#Now we have a A-Bruijn Graph with possibly cycles in it 
print("Number of glued nodes: ",glue_nodes)

#labels=dict((n,d['Base']) for n,d in G.nodes(data=True))
#nx.draw(G,labels=labels)
#plt.show()

#############################################
# Bulge Collapsing ##########################
#############################################

#First want to identify bulges in a graph



############################################
# Simplify Graph and/or find blocks ########
############################################

#Define a function to be used recursively to check for each succesor node whether it only has one in or out

def successor_check(graph,n,tmerge):
	ess = [node for node in graph.successors(n)] #Get list of succesors

        #Check for all successors
	for s in ess:
		#if(len(graph.out_edges([s])) <= 1 and len(graph.in_edges([s])) <= 1): #I.e. if only one outgoing edge and one incoming edge it belongs to same block
		if(len(graph.in_edges([s])) <= 1 and len(ess) <= 1): #Succesor node only has one incoming path and is the only option for the previous node
			tmerge.append(s)
			successor_check(graph,s,tmerge)

        #Will recursively run until there is no successor node to add then we will return the list of nodes to merge
	return(tmerge)


def merge_nodes(lnodes,graph): #Given a list of nodes merge them
        #Effectively merge all nodes onto first node

        #Redirect last node in lists edge to first node in list
	for n1,n2,d in graph.out_edges(lnodes[-1],data=True):
		graph.add_edge(lnodes[0],n2,d)

        #Get base sequence for full merge
	seq = graph.node[lnodes[0]]['Base']
	for i in range(1,len(lnodes)):
		seq = seq + graph.node[lnodes[i]]['Base']
		graph.remove_node(lnodes[i])

	#Add full sequence of merged bases to first node
	graph.node[lnodes[0]]['Base'] = seq
	return(lnodes[0]) #Return Node id of merged node


already_merged = []

#Loop through nodes

print("Simplifying Graph chains")

#Copy graph before simplifying
C = G.to_directed()
conmerge=True
if(conmerge):
	for n,d in C.nodes(data=True):

        	#if(len(C.out_edges([n])) > 1 or len(C.in_edges([n])) > 1 ): continue #If node has more than one edge coming in or out if itself skip
		if(len(C.out_edges([n])) >1 ): continue #Continue if node branches, that is to say if has more than one out edge
		if(n in already_merged): continue

		to_merge = [n]
		tmerge = successor_check(C,n,to_merge)

		if(len(tmerge) > 1):
			l = merge_nodes(tmerge,C)
			for tm in tmerge:
				already_merged.append(tm)



####################################################
####### Whirl Elimination     ######################
####################################################

whirl_removal = False
if(whirl_removal):


	#Find all whirls
	print("Finding Whirls...")
	whirls = list(nx.simple_cycles(C))
	print("DONE")

	#UNSOLVED: What if multiple interacting whirls exist?

	#Loop through each whirl
	while len(whirls) > 0:

		#Print Graph in current state
		#labels=dict((n,len(d['Base'])) for n,d in C.nodes(data=True))
		#labels=dict((n,d['Base']) for n,d in C.nodes(data=True))
		#nx.draw(C,labels=labels)
		#plt.show()

		whirl = whirls[0]
		M_node = None
		Multi = 0

		#Find Highest multiplicity node in loop to use for breaking point of cycle
		print(whirl)
		for node in whirl:
			temp = len(C.out_edges([node])) + len(C.out_edges([node]))
			if(temp >= Multi):
				Multi = temp
				M_node = node

		iM = whirl.index(M_node)
		iM = whirl[iM]

		#Make a copy of node
		C.add_node(Node_index)
		C.node[Node_index]['Base'] = C.node[iM]['Base']


		### Create edges in new node and remove some from old node
		#Copy out edges to new node and remove from old
		for n1,n2,d in C.out_edges(iM,data=True):
			if(n2 not in whirl):
				C.add_edge(Node_index,n2,d)
				C.remove_edge(iM,n2)
		

		#Get In edge to new node and remove from old
		for n1,n2,d in C.in_edges(iM,data=True):
			if(n1 in whirl): 
				C.add_edge(n1,Node_index,d)
				C.remove_edge(n1,iM)

		Node_index= Node_index+1

		#Now recalculate whirls
		whirls = list(nx.simple_cycles(C))

############################################
# Check if cyclic ##########################
############################################

if(nx.is_directed_acyclic_graph(C)):
	print("Succesfully created a DAG\n")
else:
	print("QQ\n")


###############################################
#Check we can get transcripts back from graph #
###############################################

checktran = False
if(checktran):
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
#base_order = nx.topological_sort(G)
#Will crash if there is a cycle, therefore do a try

try:
	base_order = nx.topological_sort(C)
except nx.NetworkXUnfeasible:
	print("CYCLESSSSSS!!!!!!")
	base_order = ''
seq =''

block_holder = []
#Save start and end co-ordinates of each block
coord = [0]

for index in base_order:
	seq = seq + C.node[index]['Base']
	block_holder.append(C.node[index]['Base'])
	coord.append(coord[-1] + len(C.node[index]['Base']))
print(seq)	

#Save sequence to file
superf = open('Super.fasta','w')
superf.write('>Super\n')
superf.write(seq + "\n")
superf.write("SuperBlocks: " + str(coord) + "\n")
superf.close()

#MAKE SAF file
super_saf = open('Super.saf','w')
super_saf.write('GeneID' + '\t' + 'Chr' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Strand' + '\n')
for i in range(0,len(coord)-1):
	super_saf.write('GeneX' '\t' + 'Chromo' + '\t' + str(coord[i]) + '\t' + str(coord[i+1]) + '\t' + '+' + '\n')
super_saf.close()



#Save graph to file
#nx.write_gpickle(C,"Simples.pkl")

#Save blocks to pickle file
pickle.dump(block_holder,open("SuperBlocks.pkl","wb"))

#Find the longest path in a DAG
def longest_path(G):
    dist = {} # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0]+1,v) for v in G.pred[node]] 
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node,(length,_)  = max(dist.items(), key=lambda x:x[1])
    path = []
    while length > 0:
        path.append(node)
        length,node = dist[node]
    return list(reversed(path))

#print("Longest Path ")
#lp = longest_path(G)
#seq =''
#for index in lp:
#        seq = seq + G.node[index]['Base']
#print(seq)


#################################################################
#################################################################

#labels=dict((n,len(d['Base'])) for n,d in C.nodes(data=True))
#labels=dict((n,d['Base']) for n,d in C.nodes(data=True))
#nx.draw(C,labels=labels)
#plt.show()

nx.write_gpickle(C, "test.gpickle")

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
