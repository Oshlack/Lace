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
def Split(genome):
	start_time = time.time()

	if(os.path.isfile(genome)):
		print("Creating a fasta file per gene...")

		#Parse Fasta file and split by gene
		fT = open(genome,'r')
		transcripts = {}
		geneid = {}
		transid ={}

		for line in fT:
			if(">" in line): #Name of
				tag = line.split('\n')[0].split('\r')[0].lstrip('>')
				transcripts[tag] = ''
	
				#Extract names
				trans_id = (tag.split('range')[0]).split('Gene_')[1]
				gene_id = (tag.split('range=')[1]).split(':')[0]
				geneid[tag] = gene_id
				transid[tag] =trans_id
			else:
				transcripts[tag] = transcripts[tag] + line.split('\n')[0].split('\r')[0]	

		#Make a file for each gene
		gene_list = set(geneid.values())

		for gene in gene_list:
			f = open('Genes/' + gene + '.fasta','w') 
			for tag in transcripts.keys():
				if(gene == geneid[tag]):
					#f.write('>' + transid[tag] + geneid[tag] + '\n') 
					f.write('>' +  (tag.split('refGene')[1]).split("'pad")[0] +  '\n')
					f.write(transcripts[tag]+'\n')
			f.close()

		#Now submit a BLAT job for each gene in parallel
		print("Now BLATTing each gene...")
		jobs = []

		for gene in gene_list:
			fname = 'Genes/' + gene + '.fasta'
			p = Process(target=BLATT,args=(fname,))
			jobs.append(p)
			p.start()
			#p.join()	

		# Wait for all worker processes to finish
		for p in jobs:
			p.join()	


		print("SPLIT ---- %s seconds ----" %(time.time()-start_time))

if __name__ == '__main__':

	if(len(sys.argv) != 2):
		print("Here we need an input fasta file, one for each cluster")
		sys.exit()
	else:
		genome = sys.argv[1]
		#All(genome)
		Split(genome)



