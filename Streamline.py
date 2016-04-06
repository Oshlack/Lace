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
#block_seq = [ ['T','A','C','G'], ['C','A','C','T'], ['G','G','C','T'], ['A','C'], ['C','T'], ['C']]
#tName = ['T1','T2','T3','T2','T3','T3']
#qName = ['T1','T2','T3','T1','T2','T1']
#tStart = [0,0,0,1,2,2]
#qStart = [0,0,0,1,2,2]
#transcripts = {'T1':'TACG','T2':'CACT','T3':'GGCT'}

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
Node_index=0
node_dict = {} #A place to find the index of a node
for key in transcripts:
	#node_dict[key] = np.full(len(list(transcripts[key])),-1,dtype=np.int64)
	node_dict[key] = [-1] * len(list(transcripts[key]))
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

			#Consideration 2 - Whirls from "inverted" parts
			#Check if two transcripts are disjoint (i.e. don't share any nodes) - if they are it is not problem to merge
			#shared_list = list(set(node_dict[tName[i]]).intersection(node_dict[qName[i]]))
			#If they are not, we need to check if we are making a loop
			#if(len(shared_list) > 0 ):
			#	#Now check that the base is not sequential
			#	if(qpos < node_dict[qName[i]].index(shared_list[-1])): #If the base intended to merge is sequentially before, don't merge as you can make a whirl
			#		continue

			
			

			#print("### MERGE START ###")
			#print("Started with:",qnid)
			#print("Merging Node: ", G.node[qnid])
			#print("To: ", G.node[tnid])
			#print("Which is:",tnid)

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

			#Recursive check that no element in node dict contains the old node which is removed
			for key in node_dict:
				if(key ==  qName[i] or key == tName[i]): continue
				if(qnid in node_dict[key]):
					node_dict[key][node_dict[key].index(qnid)] = tnid


################################
# Check whirl potential ########
################################
whirl_flag = False
#for key in transcripts:
#        if(len(node_dict[key]) != len(set(node_dict[key]))):
#                whirl_flag = True

if(whirl_flag):
        print("Whirl detected....exiting")
        sys.exit()




###################################################
## Add Edges between adjacent transcript nodes ####
###################################################
for key in node_dict:
        for j in range(0,len(node_dict[key])-1):
                G.add_edge(node_dict[key][j],node_dict[key][j+1])


############################################
# Simplify Graph and/or find blocks ########
############################################

#Define a function to be used recursively to check for each succesor node whether it only has one in or out

def successor_check(G,n,tmerge):
        ess = [node for node in G.successors(n)]

        #Check for all successors
        for s in ess:
                if(len(G.out_edges([s])) == 1 and len(G.in_edges([s])) == 1): #I.e. if only one outgoing edge and one incoming edge it belongs to same block
                        tmerge.append(s)
                        successor_check(G,s,tmerge)

        #Will recursively run until there is no successor node to add then we will return the list of nodes to merge
        return(tmerge)


def merge_nodes(lnodes,graph): #Given a list of nodes merge them
	print(lnodes)
        #Effectively merge all nodes onto first node

        #Redirect last node in lists edge to first node in list
	for n1,n2 in graph.out_edges(lnodes[-1]):
		graph.add_edge(lnodes[0],n2)

        #Get base sequence for full merge
	seq = graph.node[lnodes[0]]['Base']
	for i in range(1,len(lnodes)):
		seq = seq + graph.node[lnodes[i]]['Base']
		graph.remove_node(lnodes[i])

	#Add full sequence of merged bases to first node
	graph.node[lnodes[0]]['Base'] = seq
	return(lnodes[0]) #Return Node id of merged node


#already_merged = []

#Loop through nodes

#Copy graph before simplifying
C = G.to_directed()

#for n,d in G.nodes(data=True):

#        if(len(G.out_edges([n])) > 1 or len(G.in_edges([n])) > 1 ): continue #If node has more than one edge coming in or out if itself skip
#        if(n in already_merged): continue

#        to_merge = [n]
#        tmerge = successor_check(G,n,to_merge)


#        if(len(tmerge) > 1):
#                l = merge_nodes(tmerge)
#                for tm in tmerge:
#                        already_merged.append(tm)



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
	for whirl in whirls:
		M_node = None
		Multi = 0

		#Find Highest multiplicity node in loop to use for breaking point 
		for node in whirl:
			temp = len(C.out_edges([node])) + len(C.out_edges([node]))
			if(temp > Multi):
				Multi = temp
				M_node = node

		
		print(whirl)
		iM = whirl.index(M_node)
		print(iM)

		#Merge all nodes, after highest multiplicity node together
		if(len(whirl[(iM+1):]) == 1): wnode = M_node
		elif(len(whirl[(iM+1):]) == 0): 
			wnode = M_node
			ppnode = merge_nodes(whirl[0:iM],C)

		else:
			wnode = merge_nodes(whirl[(iM+1):],C)

			#Merge all nodes in whirl up to (and including)
			ppnode = merge_nodes(whirl[0:iM+1],C)

		#Make a copy of ppnode
		C.add_node(Node_index)
		C.node[Node_index]['Base'] = C.node[ppnode]['Base']


		### Create edges in new node and remove some from old node
		#Copy out edges to new node and remove from old
		for n1,n2 in C.out_edges(ppnode):
			if(n2 not in whirl):
				C.add_edge(Node_index,n2)
				C.remove_edge(ppnode,n2)
		

		#Get In edge to new node and remove from old
		for n1,n2 in C.in_edges(ppnode):
			if(n1 in whirl): 
				C.add_edge(n1,Node_index)
				C.remove_edge(n1,ppnode)

		#Now add edge from break to whirl
		C.add_edge(ppnode,node)
		C.add_edge(wnode,Node_index)

		
			
		Node_index= Node_index+1

###############################################
#Check we can get transcripts back from graph #
###############################################

#for key in transcripts:
#	path_alternatives = []
	#print("Transcript: ", transcripts[key])
#	pathl = nx.all_simple_paths(G,node_dict[key][0],node_dict[key][-1])
#	for path in pathl:
#		seq= ''
#		#print(path)
#		for b in path:
#			#print(G.node[b]['Base'])
#			seq = seq + G.node[b]['Base']
#		#print(seq)
#		path_alternatives.append(seq)

#	#Loop through and check the strings match
#	for path in path_alternatives:
#		if(path == transcripts[key]): print("Found path for: ", key)
		


#Order base in graph
#base_order = nx.topological_sort(G)
#Will crash if there is a cycle, therefore do a try
try:
	base_order = nx.topological_sort(G)
except nx.NetworkXUnfeasible:
	print("CYCLESSSSSS!!!!!!")
	sys.exit()

seq =''
for index in base_order:
	seq = seq + G.node[index]['Base']
print(seq)	

#Save sequence to file
superf = open('Super.fasta','w')
superf.write('>Super\n')
superf.write(seq)
superf.close()


#Save graph to file
#nx.write_gpickle(C,"Simples.pkl")


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

#labels=dict((n,d['Base']) for n,d in C.nodes(data=True))
#nx.draw(C,labels=labels)
#plt.show()

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
