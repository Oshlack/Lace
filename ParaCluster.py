#The idea here is to build a script which takes a list of fasta files (each file being a set of transcripts from a Corset cluster which supposedly come from the same gene).
#Then given each cluster, there should be a seperate parallelised job for each file which constructs the SuperTranscript for that cluster.

import multiprocessing, logging
from multiprocessing import Pool
from multiprocessing import Process
import os
from BuildSuperTranscript import SuperTran
import sys
import time

#logger = multiprocessing.log_to_stderr()
#logger.setLevel(multiprocessing.SUBDEBUG)

def info(title):
	print(title)
	print('module name:',__name__)
	if hasattr(os,'getppid'): #only avaliable on Unix
		print('parent process:',os.getppid())
	print('process id:',os.getpid())

def worker(fname):
	seq =''
	ann = ''
	try:
		seq,ann = SuperTran(fname)
	except:
		print("Failed:", fname)
	return seq,ann

def Para(clustlist):
	clusters = []
	if(os.path.isfile(clustlist)):
		f = open(clustlist,'r')
		for line in f:
			clusters.append(line.split("\n")[0])

	else:
		print("Input list not found, exiting...")
		sys.exit()

	# BY POOL
	ncore = multiprocessing.cpu_count()
	pool = Pool(processes=ncore)
	result= pool.map_async(SuperTran,clusters) 
	pool.close()

	results = result.get()

	write_dir = os.path.dirname(clustlist)
	if(write_dir == ''): write_dir = '.'

	#Write Overall Super Duper Tran
	superf = open(write_dir + '/' +'SuperDuper.fasta','w')
	supgff = open(write_dir + '/' +'SuperDuper.gff','w')

	for i,clust in enumerate(clusters):
		superf.write('>' + clust + '\n')
		superf.write(result[i][0] + '\n')

	#Write Super gff
	for res in result:
        	supgff.write(res[1])
		

#Split fasta file into genes first then parallelise the BLAT for the different genes
def Split(genome,corsetfile):
	start_time = time.time()

	#Find working directory
	dir = os.path.dirname(corsetfile)
	if(dir==''): dir='.'


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
			clust = clust.replace('/','-')
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
		if('None' in gene_list): gene_list.remove('None') #Remove the placer holder for the ' None' which were transcripts not mapped to clusters in corset

		cnts = []
		for gene in gene_list:

			#Count number of transcripts assigned to a cluster
			cnt =0
			for val in cluster.values():
				if(val == gene): cnt = cnt + 1
			cnts.append(cnt)

			fn = dir + '/' + gene + '.fasta' #General		

			if(os.path.isfile(fn)): continue	#If already file
			

			f = open(fn,'w') 
			for tag in transcripts.keys():
				if(gene == geneid[tag]):
					f.write('>' + tag +  '\n')
					f.write(transcripts[tag]+'\n')
			f.close()

		#Now submit Build Super Transcript for each gene in parallel
		print("Now BLATTing each gene...")
		jobs = []

		fnames = []
		for gene in gene_list:
			fname = dir + '/' + gene + '.fasta'
			fnames.append(fname)

		# BY POOL
		ncore = multiprocessing.cpu_count()
		pool = Pool(processes=ncore)
		results = pool.map_async(worker,fnames)
		pool.close()
		pool.join()
		result = results.get()


		#Write Overall Super Duper Tran
		superf = open(dir + '/' +'SuperDuper.fasta','w')
		supgff = open(dir + '/' +'SuperDuper.gff','w')

		for i,clust in enumerate(fnames):
			superf.write('>' + clust + ' Number of transcripts: ' + str(cnts[i]) +  '\n')
			superf.write(result[i][0] + '\n')

		#Write Super gff
		for res in result:
                        supgff.write(res[1])

		print("SPLIT ---- %s seconds ----" %(time.time()-start_time))
	
if __name__ == '__main__':

	if(len(sys.argv) == 1):
		print("Here we need an input fasta file and an input corset file")
		print("Or a list of fasta files in a text file")
		sys.exit()

	else:
		if(len(sys.argv) == 3):
			genome = sys.argv[1]
			corset = sys.argv[2]
			Split(genome,corset)

		if(len(sys.argv) == 2):
			clustlist = sys.argv[1]
			Para(clustlist)
			#Ser(clustlist)
