# Author: Anthony Hawkins

#A script to construct a SuperTranscript given a .fasta file of transcript sequence
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
import traceback

sys.setrecursionlimit(10000)


#Define a function to be used recursively to check for each succesor node whether it only has one in or out
def successor_check(graph,n,tmerge):
    ess = [node for node in graph.successors(n)] #Get list of succesors

    #Check for all successors
    for s in ess:
        if(len(graph.in_edges([s])) <= 1 and len(ess) <= 1): #Succesor node only has one incoming path and is the only option for the previous node
            tmerge.append(s)
            successor_check(graph,s,tmerge)

    #Will recursively run until there is no successor node to add then we will return the list of nodes to merge
    return(tmerge)


def merge_nodes(lnodes,graph): #Given a list of nodes merge them
    #Effectively merge all nodes onto first node

    #Redirect last node in lists edge to first node in list
    for out_edge in graph.out_edges(lnodes[-1],data=True):
        n1,n2,d = out_edge
        if d:
            graph.add_edge(lnodes[0],n2,d)
        else:
            # d looks empty, so removing it to avoid: TypeError: add_edge() takes exactly 3 arguments (4 given)
            graph.add_edge(lnodes[0],n2)
    
    #Get base sequence for full merge
    seq = graph.node[lnodes[0]]['Base']
    for i in range(1,len(lnodes)):
        seq = seq + graph.node[lnodes[i]]['Base']
        graph.remove_node(lnodes[i])

    #Add full sequence of merged bases to first node
    graph.node[lnodes[0]]['Base'] = seq
    return(lnodes[0]) #Return Node id of merged node

#Given a transcript find its reverse compliment
def Reverse_complement(transcript):
    rev_tran = []
    for base in transcript:
        rev_base = ''
        if(base.lower()=='a'): rev_base = 't'
        elif(base.lower()=='c'): rev_base = 'g'
        elif(base.lower()=='t'): rev_base = 'a'
        else: rev_base = 'c'
        rev_tran.append(rev_base)
    return ('').join(rev_tran[::-1]) #Now return the list reversed


#Define direction of transcripts in fasta file
#Loop through table from psl
#Remove any repeated rows from table and define directionality of transcripts
def filt_dir(table):
    pair_list = [] # Which pairs of transcripts have been blatted
    rem_row = [] #A list of rows to remove from dataframe
    trandir = {} #A dictionary holding the directionality of transcripts, this is defined arbitrarily with respect to whatever comes first in psl file


    for i in range(0,len(table)):

        ########################
        # Filtering ############
        ########################

        #Don't allign the transcripts against each other twice...
        #I.e. BLAT does T1 vs T2 but also T2 vs T1 (which should be the same give or take)

                #TName + QName
        paired = table.iloc[i,13] + table.iloc[i,9]
        #If that pair of transcripts is already in the table then remove from table
        if(paired in pair_list):
            rem_row.append(i)
            continue
        else: pair_list.append(table.iloc[i,9]+table.iloc[i,13])

        #Obviously we do not need the rows where we blat a transcript against itself
        if(table.iloc[i,13] == table.iloc[i,9]):
            rem_row.append(i)
            continue

        ##################
        # Directionality #
        ##################

        #Directionality of transcripts
        tName = table.iloc[i,13]
        qName = table.iloc[i,9]
        strand = table.iloc[i,8]

        #Check if direction of one of the transcripts is defined
        if((tName in trandir) and (qName in trandir)): #Directionality of both transcripts already defined
            continue

        #Transcript directiionality already defined
        elif(tName in trandir): 
            tdir = trandir[tName]
            if(tdir == strand): #If both the transcript and strand same
                trandir[qName] = '+'
            else:
                trandir[qName] = '-'

        #Query transcript already defined
        elif(qName in trandir): 
            qdir = trandir[qName]
            if(qdir == strand): #If both the query and strand same
                                trandir[tName] = '+'
            else:
                trandir[tName] = '-'
        
        #Neither of the sequences defined
        else:
            if(strand=='+'): #Both transcripts in the same direction
                trandir[qName] = '+'
                trandir[tName] = '+'

            else: #arbitrarily define transcript as postive and query as negative
                trandir[tName] = '+'
                trandir[qName] = '-'
        
    ###############
    # Removal #####
    ###############
    #Remove rows where the transcript pair is repeated
    table = table.drop(table.index[rem_row])        

    return(table,trandir)
                    

    

#Main function to produce a SuperTranscript
def SuperTran(fname,verbose=True):
    
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

    transcript_status = len(transcripts)
    whirl_status = 0

    #If there is only one transcript in this file, then simply that transcript is the super transcript...
    if(len(transcripts) == 1):
        if(verbose): print("One\n") 
        seq = next(iter(transcripts.values())) #Python 3 specific codee...
        anno = (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + '1' + '\t' + str(len(seq)) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + '.'  + '\n'

    else:
        #Try topo sorting a graph
        try:
            seq, anno, whirl_status  = BuildGraph(fname,transcripts,verbose)

        except: #Graph building failed, just take longest transcript (or concatenate all transcripts)
            traceback.print_exc()
            temp = 0
            seq = ''
            print('FAILED to construct')
            for val in list(transcripts.values()):
                if(len(val) > temp):
                    temp = len(val)
                    seq = ''.join(val)

            anno = (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + '1' + '\t' + str(len(seq)) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' + '.' + '\n'
    
    return(seq,anno,whirl_status,transcript_status)

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
        BLAT_command = "blat %s %s -minIdentity=98  %s.psl" %(fname,fname,fname.split('.fasta')[0]) #This gets almost exact matches
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
    trandir = {} # A dictionary defining the directionality of transcripts

    #Filter psl table where two transcripts can have multiple rows (usually because blat does T1 vs T2 then T2 vs T1 later)
    bData, trandir = filt_dir(bData)

    if(len(bData['strand'].unique()) > 1 or bData['strand'].unique() == '-'): #That is we have both pos and neg strands or potentially two transcripts with a negative strand
        print("Double Stranded Contigs\n")

        #Re-correct the transcripts to be the reverse compliments if one of the transcripts has a negative directionality
        for key in trandir:
            if(trandir[key] == '-'):
                transcripts[key] = Reverse_complement(transcripts[key])

    
        #Write a fasta file with the contigs all in same direction
        fcorr = fname.split('.fasta')[0] + "_stranded" + ".fasta"
        print(fcorr)
        fc = open(fcorr,"w")
        for key in transcripts:
            fc.write(">" + key + "\n")
            fc.write(transcripts[key]+"\n")
        fc.close()    

        #Re-BLAT
        reblat = "blat %s %s -maxGap=0 -minIdentity=98  %s.psl" %(fcorr,fcorr,fcorr.split('.fasta')[0]) #This gets almost exact matches
        os.system(reblat)
        bData = pd.read_table(fcorr.split('.fasta')[0] + '.psl',sep='\t',header = None,names=Header_names,skiprows=5)
    
        #Re-filter
        bData, trandir = filt_dir(bData)

    
    for i in range(0,len(bData)):
        
        #Extract the info
        seq=list(transcripts[bData.iloc[i,9]]) #Get sequence from query name
        block_sizes = (bData.iloc[i,18]).rstrip(',').split(',')
        qStarts = (bData.iloc[i,19]).rstrip(',').split(',')
        tStarts = (bData.iloc[i,20]).rstrip(',').split(',')

        for j in range(0,len(qStarts)):
            block_seq.append(seq[int(qStarts[j]):(int(qStarts[j])+int(block_sizes[j]))]) #This is purely used for the size of the sequence
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
    reverse_node_dict = []
    for key in transcripts:
        node_dict[key] = [-1] * len(list(transcripts[key]))

    #############################
    # ADD NODES AND EDGES #######
    #############################

    if(verbose): print("Setting up graph: adding a node for each base")

    #Add a node for every base in each transcript
    for key in transcripts:
        for i in range(0,len(transcripts[key])):
            #Add node to graph
            G.add_node(Node_index)

            #Add attributes
            G.node[Node_index]['Base'] = transcripts[key][i]

            #Add to the reverse dictionary
            reverse_node_dict.append( {key : [i]} )

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

    if(verbose): print("Add Edges between adjacent transcript nodes")

    for key in node_dict:
        for j in range(0,len(node_dict[key])-1):
            G.add_edge(node_dict[key][j],node_dict[key][j+1])


    ####################################################
    ## Let the gluing commence #########################
    ####################################################
    if(verbose): 
        print("Merging nodes")
        
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
                if(qName[i] in reverse_node_dict[tnid]):
                    if(reverse_node_dict[tnid][qName[i]]!=[qpos]): continue

                #If the node you are intending to merge is already merged to somewhere else on the transcript string, dont merge as can cause wirls
                if(tName[i] in reverse_node_dict[qnid]): continue

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

                #Remove query node since we have now merged it to the transcript node
                G.remove_node(qnid)

                #Change Dictionary Call for query node and 
                #Check that no element in node dict contains the old node which is removed, if it does replace it
                for key in reverse_node_dict[qnid]:
                    for pos in reverse_node_dict[qnid][key]:
                        node_dict[key][pos] = tnid

                #Update the reverse dictionary
                for key in reverse_node_dict[qnid]:
                    if key in reverse_node_dict[tnid]:
                        reverse_node_dict[tnid][key].extend(reverse_node_dict[qnid][key])
                        reverse_node_dict[tnid][key]=set(reverse_node_dict[tnid][key])
                    else:
                        reverse_node_dict[tnid][key]=reverse_node_dict[qnid][key]
                reverse_node_dict[qnid]={}


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
        c_nodes = C.nodes(data=True)
        # avoid: RuntimeError: dictionary changed size during iteration w/ py3.6.5
        ns = list()
        for n,d in c_nodes:
            ns.append(n)

        for n in ns:
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
    whirl_status = 0
    if(whirl_removal):
        #Find all whirls
        #print("Finding Whirls...")
        whirls = list(nx.simple_cycles(C))
        #print("DONE")
        whirl_status = len(whirls) #Report initial numbe of whirls in graph

        #Loop through each whirl
        while len(whirls) > 0:
            whirl = whirls[0]
            M_node = None
            Multi = 100000000

                        #Find the node with smallest sequence and break there (instead of the highest multiplicty)
            for node in whirl:
                temp = len(C.node[node]['Base'])
                if(temp <= Multi):
                    Multi =temp
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
        if(verbose): print("CYCLESSSSSS!!!!!!")
        anno = ''
        return("CYCLE",anno)
    
        
    seq =''
    coord = [0]
    for index in base_order:
        seq = seq + C.node[index]['Base']
        coord.append(coord[-1] + len(C.node[index]['Base'])) #0-based co-ordinates

    #String for annotation file
    anno = ''
    for i in range(0,len(coord)-1):
        anno = anno + (fname.split('/')[-1]).split('.fasta')[0] + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(coord[i]+1) + '\t' + str(coord[i+1]) + '\t' + '.' + '\t' +'.' + '\t' + '0' + '\t' +  '.'  + '\n' #GFF2 format - 1 base for igv

    return(seq,anno,whirl_status)

if __name__ == '__main__':
    ''' Takes one fasta file which contains all transcripts in cluster (gene) and builds a super transcript from it, outputing the sequence'''

    if(len(sys.argv) != 2):
        print('Function takes one fasta file as input')
        exit

    else:
        fname = sys.argv[1]
        seq,anno,whirl_status,transcript_status = SuperTran(fname,verbose=True)

        print(seq)
        print(anno)
        print(whirl_status)
        print(transcript_status)

    
