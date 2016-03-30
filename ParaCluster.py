#The idea here is to build a script which takes a list of fasta files (each file being a set of transcripts from a Corset cluster which supposedly come from the same gene).
#Then given each cluster, there should be a seperate parallelised job for each file which constructs the SuperTranscript for that cluster.

from multiprocessing import Pool
from multiprocessing import Process
import os
from BuildSuperTranscript import SuperTran
import sys

def info(title):
	print(title)
	print('module name:',__name__)
	if hasattr(os,'getppid'): #only avaliable on Unix
		print('parent process:',os.getppid())
	print('process id:',os.getpid())

def f(name):
	info('function f')
	print('hello',name)


def Para(clustlist):
	clusters = []
	if(os.path.isfile(clustlist)):
		f = open(clustlist,'r')
		for line in f:
			clusters.append(line.split("\n")[0])

	else:
		print("Input list not found, exiting...")
		sys.exit()

	jobs = []
	for cluster in clusters:
		p = Process(target=SuperTran,args=(cluster,))
		jobs.append(p)
		p.start()
		p.join()


if __name__ == '__main__':

	if(len(sys.argv) != 2):
		print("Here we need an input text file with the list of all fasta files, one for each cluster")
		sys.exit()

	else:
		clustlist = sys.argv[1]
		Para(clustlist)
