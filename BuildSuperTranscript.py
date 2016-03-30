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

sys.setrecursionlimit(10000) # 10000 is an example, try with different values


def SuperTran(fname):

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


	###################################################
	# Loop pairwise through transcripts and BLAT allign
	###################################################

	#This may well be changed/ skipped/ parallelised
	if(os.path.isfile('outall.psl') == False):
		BLAT_command = "./blat %s %s -maxGap=0 -minIdentity=100 -maxIntron=0 outall.psl" %(fname,fname) #This gets almost exact matches
		os.system(BLAT_command)


	#First read in BLAT output:
	Header_names = ['match','mismatch','rep.','N\'s','Q gap count','Q gap bases','T gap count','T gap bases','strand','Q name','Q size','Q start','Q end','T name','T size','T start','T end','block count','blocksizes','qStarts','tStarts']

	bData = pd.read_table('ESR1_exakt.psl',sep='\t',header = None,names=Header_names,skiprows=5)

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


	########################
	# Construct the Graph ##
	########################
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


			#Co-ordinate for base in transcript and query sequences
			tpos = j + tStart[i]
			qpos = j + qStart[i]
		
			#Node index for q base and t base
			qnid = node_dict[qName[i]][qpos] #Node index in graph for base in query
			tnid = node_dict[tName[i]][tpos] #Node index in graph for base in transcript
			

			#Check if node already there either with the transcript id or query id

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




	#Order base in graph
	base_order = nx.topological_sort_recursive(G)
	seq =''
	for index in base_order:
		seq = seq + G.node[index]['Base']

	#Save sequence to file
	superf = open('Super.fasta','w')
	superf.write('>Super\n')
	superf.write(seq)
	superf.close()

	print("---- %s seconds ----" %(time.time()-start_time))	
	print(seq)	
	
	return(seq)

if __name__ == '__main__':
	''' Takes one fasta file which contains all transcripts in cluster (gene) and builds a super transcript from it, outputing the sequence'''
	

	if(len(sys.argv) != 2):
		print('Function takes one fasta file as input')
		exit

	else:
		fname = sys.argv[1]
		seq = SuperTran(fname)
		

	
