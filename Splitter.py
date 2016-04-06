#A script for taking in a fasta file and splitting it up by gene id and then running BLAT on them, or running BLAT on the large fasta file and then afterwards grouping.

import time
import sys
import os
from multiprocessing import Pool
from multiprocessing import Process


#Call BLAT just once for every transcript against every other transcript
def All(genome):
	start_time = time.time()
	if(os.path.isfile(genome)):
		BLAT_command = "./blat %s %s -maxGap=0 -minIdentity=100 -maxIntron=0 gene.psl" %(genome,genome) #This gets almost exact matches
		os.system(BLAT_command)

	print("ALL ---- %s seconds ----" %(time.time()-start_time))


#Pairwise allign all transcripts in file
def BLATT(fname):
	if(os.path.isfile(fname)):
		BLAT_command = "./blat %s %s -maxGap=0 -minIdentity=100 -maxIntron=0 %s.psl" %(fname,fname,fname.split('.fasta')[0]) #This gets almost exact matches
		os.system(BLAT_command)	
	
	else:
		print("Input file not found, exiting...")
		sys.exit()





#Split fasta file into genes first then parallelise the BLAT for the different genes
def Split(genome,corsetfile):
	start_time = time.time()


	#First create dictionary from corset output assigning a transcript to a cluster
	#Note: Corset throws away transcripts which don't have many read (> 10) to them
	#So we need to also ignore them

	cluster = {}
	if(os.path.isfile(corsetfile)):
		print("Creating dictionary of transcripts in clusters...")
		corse = open(corsetfile,'r')
		for line in corse:
			tran = line.split()[0]	
			clust = line.split()[1].rstrip('/n')
			cluster[tran] = clust			


	#Now loop through fasta file
	if(os.path.isfile(genome)):
		print("Creating a fasta file per gene...")
		#We only want to include transcripts deemed worthy by Corset

		#Parse Fasta file and split by gene
		fT = open(genome,'r')
		transcripts = {}
		geneid = {}
		transid ={}

		for line in fT:
			if(">" in line): #Name of
				tag = (line.split()[0]).lstrip('>')

				transcripts[tag] = ''
	
				#Assign names
				if(tag in cluster.keys()): geneid[tag] = cluster[tag] #If assigned by corset
				else: geneid[tag] = 'None'
				transid[tag] = tag
			else:
				transcripts[tag] = transcripts[tag] + line.split('\n')[0].split('\r')[0]	

		#Make a file for each gene
		gene_list = set(geneid.values())
		gene_list.remove('None') #Remove the placer holder for the ' None' which were transcripts not mapped to clusters in corset

		for gene in gene_list:

			#Count number of transcripts assigned to a cluster
			#cnt =0
			#for val in cluster.values():
			#	if(val == gene): cnt = cnt + 1

			f = open(corsetfile.split('clusters')[0] + gene + '.fasta','w') 
			for tag in transcripts.keys():
				if(gene == geneid[tag]):
					f.write('>' + tag +  '\n')
					f.write(transcripts[tag]+'\n')
			f.close()

		#Now submit a BLAT job for each gene in parallel
		print("Now BLATTing each gene...")
		jobs = []

		fnames = []
		for gene in gene_list:
			fname = corsetfile.split('clusters')[0] + gene + '.fasta'
			fnames.append(fname)

		# BY POOL
		pool = Pool(processes=4)
		result = pool.map(BLATT,fnames)
		pool.close()
		pool.join()

		

		print("SPLIT ---- %s seconds ----" %(time.time()-start_time))

if __name__ == '__main__':

	if(len(sys.argv) != 3):
		print("Here we need an input fasta file and an input corset file")
		sys.exit()
	else:
		genome = sys.argv[1]
		clusters = sys.argv[2]
		#All(genome)
		Split(genome,clusters)



