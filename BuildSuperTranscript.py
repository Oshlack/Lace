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
import os
from matplotlib.pyplot import cm 

sys.setrecursionlimit(10000)


	#Define a function to be used recursively to check for each succesor node whether it only has one in or out

def successor_check(graph,n,tmerge):
	ess = [node for node in graph.successors(n)] #Get list of succesors

	#Check for all successors
	for s in ess:
		#if(len(graph.out_edges([s])) == 1 and len(graph.in_edges([s])) == 1): #I.e. if only one outgoing edge and one incoming edge it belongs to same block
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


def SuperTran(fname,verbose=False):
	
	#Start Clock for timing
	start_time = time.time()

	#####################
	#Read in transcripts
	#####################
	if(not os.path.isfile(fname)):
		print("File name corrupt")
		sys.exit()

	fT = open(fname,'r')
	transcripts = {}
	tName = ''

	for line in fT:
		if(">" in line): #Name of 
			tName = line.split('\n')[0].split('\r')[0].lstrip('>')
			transcripts[tName] = ''
		
		else:
			transcripts[tName] = transcripts[tName] + line.split('\n')[0].split('\r')[0]

	seq = ''

	#If there is only one transcript in this file, then simply that transcript is the super transcript
	if(len(transcripts) == 1):
		if(verbose): print("One\n") 
		seq = next(iter(transcripts.values())) #Python 3 specific codee...
		#anno = (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'Chromo' + '\t' + '0' + '\t' + str(len(seq)) + '\t' + '+' + '\n'
		anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + 0 + '\t' + str(len(seq)) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + '.'  + '\n'

	else:
		#Try topo sorting a graph
		try:
			seq, anno = BuildGraph(fname,transcripts,verbose)

		except: #Graph building failed, just take longest transcript or (concatenate all transcripts)
			temp = 0
			seq = ''
			print('FAILED to build graph and topo sort')
			for val in transcripts:
				if(len(val) > temp):
					temp = len(val)
					seq = ''.join(val)

			#anno = (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'Chromo' + '\t' + '0' + '\t' + str(len(seq)) + '\t' + '+' + '\n'
			anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + 0 + '\t' + str(len(seq)) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + '.' + '\n'
	
	#print("---- %s seconds ----" %(time.time()-start_time))
	return(seq,anno)

def BuildGraph(fname,transcripts,verbose=False):
	# A Function to build a bruijn graph/splice node graph based on the transcripts assigned to a given cluster
	# One node per base in the transcript are created, then based on pairwise allignments of transcripts (using BLAT)
	# nodes in overlapping transcripts are glued together
	# Then the graph is collapsed into superblocks where each node is built of a collapsed chain of nodes with one incoming and outgoing edge
	# Finally a topological sorting is made 

	###################################################
	# Loop pairwise through transcripts and BLAT allign
	###################################################

	#This may well be changed/ skipped/ parallelised
	if(not os.path.isfile(fname.split('.fasta')[0] + '.psl')):
		BLAT_command = "./blat %s %s -maxGap=0 -minIdentity=100 -maxIntron=0 %s.psl" %(fname,fname,fname.split('.fasta')[0]) #This gets almost exact matches
		os.system(BLAT_command)


	#First read in BLAT output:
	Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']

	bData = pd.read_table(fname.split('.fasta')[0] + '.psl',sep='\t',header = None,names=Header_names,skiprows=5)

	###############
	#Now extract the sequences from the blocks using the coordinates from BLAT
	###############

	block_seq= []
	tName = [] #The name of the transcript to which the block is shared in 
	qName = []
	tStart = [] #Start co-ordinate of the block in the transcript
	qStart = [] #Start co-ordinate of the block in query
	pair_list = [] # Which pairs of transcripts have been blatted

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
		
		#Extract thhe info
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

	if(verbose): print("Constructing and merging nodes in graphs based on Blocks and Transcripts")


	########################
	# Construct the Graph ##
	########################
	G= nx.DiGraph()
	Node_index=0
	node_dict = {} #A place to find the index of a node
	for key in transcripts:
		node_dict[key] = [-1] * len(list(transcripts[key]))

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

	###################################################
        ## Add Edges between adjacent transcript nodes ####
        ###################################################
	for key in node_dict:
		for j in range(0,len(node_dict[key])-1):
			G.add_edge(node_dict[key][j],node_dict[key][j+1])


	####################################################
	## Let the gluing commence #########################
	####################################################

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

			#If they are not the same node, we need to merge them and add the same edges, redirect the query node to the transcript node
			if(qnid != tnid): 

				#Consideration - Whirls from repeated sections
				#Check if transcript node id already used for another base on the query string
				try:
					ll = node_dict[qName[i]].index(tnid)

				except ValueError:
					ll = -1

				#Check whether the transcript node you are merging to, isnt already in the query string
				if(ll >= 0 and ll != qpos): continue

				#If the node you are intending to merge is already merged to somewhere else on the transcript string, dont merge as can cause wirls
				if(qnid in node_dict[tName[i]]): continue			


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

				
				#Recursive check that no element in node dict contains the old node which is removed, if it does replace it...perhaps think of another way...
				for key in node_dict:
					if(qnid in node_dict[key]):
						node_dict[key][node_dict[key].index(qnid)] = tnid


	if(verbose): print("Adding in node edges based on transcripts")


	############################################
	# Simplify Graph and/or find blocks ########
	############################################

	already_merged = []

	#Loop through nodes
	if(verbose): print("Simplifying Graph chains")

	#Copy graph before simplifying
	C = G.to_directed()
	conmerge=True
	if(conmerge == True):
		for n,d in C.nodes(data=True):
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
	whirl_removal = True
	if(whirl_removal):
		#Find all whirls
		#print("Finding Whirls...")
		whirls = list(nx.simple_cycles(C))
		#print("DONE")


		#Loop through each whirl
		while len(whirls) > 0:

			whirl = whirls[0]
			M_node = None
			Multi = 0

			#Find Highest multiplicity node in loop to use for breaking point of cycle
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

	

	#Will crash if there is a cycle, therefore do a try
	try:
		base_order = nx.topological_sort(C)

	except nx.NetworkXUnfeasible:
		print("CYCLESSSSSS!!!!!!")
		anno = ''
		return("CYCLE",anno)
	
		
	seq =''
	coord = [0]
	for index in base_order:
		seq = seq + C.node[index]['Base']
		coord.append(coord[-1] + len(C.node[index]['Base']))

	#String for annotation file
	anno = ''
	for i in range(0,len(coord)-1):
		#anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + str(coord[i]) + '\t' + str(coord[i+1]) + '\t' + '+' + '\n' #SAF format
		#anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + str(coord[i]) + '\t' + str(coord[i+1]) + '\n' #Basic BED format
		anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(coord[i]) + '\t' + str(coord[i+1]) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' +  '.'  + '\n' #GFF2 format


	#Save sequence to file
	#superf = open('Super.fasta','w')
	#superf.write('>Super' + '\n')
	#superf.write(seq)
	#superf.close()

	return(seq,anno)

if __name__ == '__main__':
	''' Takes one fasta file which contains all transcripts in cluster (gene) and builds a super transcript from it, outputing the sequence'''
	

	if(len(sys.argv) != 2):
		print('Function takes one fasta file as input')
		exit

	else:
		fname = sys.argv[1]
		seq,anno = SuperTran(fname,verbose=True)

		print(seq)
		print(anno)
		

	
