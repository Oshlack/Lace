#!/usr/bin/env python
#Author: Anthony Hawkins

#Lacing together different transcripts into a SuperTranscript
#This is the main script which parrallelises making a SuperTranscript for each gene/cluster
#The main inputs are .fasta file containing all transcripts in all the genes/cluster you wish to constuct
#and a tab delimited text file of two columns with the mapping of transcripts <-> gene


import multiprocessing, logging
from multiprocessing import Pool
from multiprocessing import Process
import os
from Lace.BuildSuperTranscript import SuperTran
from Lace.BuildSuperTranscript import get_annotation_line
import sys
import time
import argparse
from Lace.Checker import Checker
import traceback

def worker(fname):
    print("Processing cluster %s - %s" %(fname[1],fname[0]))
    seq =''
    ann = ''
    whirl_status=0
    transcript_status=0
    try:
        seq,ann,whirl_status,transcript_status = SuperTran(fname[0])
    except:
        traceback.print_exc()
        print("Failed:", fname[0])
    return seq,ann,whirl_status,transcript_status

#A little function to move all .fasta and .psl files created into a sub directory to tidy the space
def Clean(clusters,outdir,full_clean):
    
    print("Cleaning up")

    if(not full_clean):
        mcom_mkdir = 'mkdir %s/SuperFiles' %(outdir)
        os.system(mcom_mkdir)
    for clust in clusters:
        if(full_clean):
            mcom = 'rm %s/%s.fasta %s/%s_stranded.fasta %s/%s.psl %s/%s_stranded.psl 2>/dev/null' \
                %(outdir,clust,outdir,clust,outdir,clust,outdir,clust)
        else:
            mcom = 'mv %s/%s.fasta %s/%s_stranded.fasta %s/%s.psl %s/%s_stranded.psl %s/SuperFiles 2>/dev/null' \
                %(outdir,clust,outdir,clust,outdir,clust,outdir,clust,outdir)
        os.system(mcom)


#Split fasta file into genes first then parallelise the BLAT for the different genes
def Split(genome,corsetfile,ncore,maxTran,outdir,full_clean):
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
            cluster[tran] = clust            

    #Check which clusters only contain one transcript
    single_cluster=[]
    for clus in set(cluster.values()) :
        trans=list(cluster.values()).count(clus)
        if trans==1 :
            single_cluster.append(clus)

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
                if(tag in cluster.keys() and cluster[tag] not in single_cluster ): #remove single transcript clusters
                    geneid[tag] = cluster[tag] #If assigned by corset
                else: geneid[tag] = 'None'
                transid[tag] = tag
            else:
                transcripts[tag] = transcripts[tag] + line.split('\n')[0].split('\r')[0]    

        #Make a file for each gene
        gene_list = set(geneid.values())
        if('None' in gene_list): gene_list.remove('None') #Remove the placer holder for the ' None' which were transcripts not mapped to clusters in corset

        cnts = []
        for gene in gene_list:
            #Count number of transcripts assigned to cluster
            cnt=0
            for val in cluster.values():
                if(val ==gene):cnt=cnt+1
            cnts.append(cnt)
            if(cnt > maxTran): print("WARNING: Lace will only take the first " + str(maxTran) +" transcripts since there are too many transcripts in "+gene) 
        
            fn = outdir + '/' + gene + '.fasta' #General        
            if(os.path.isfile(fn)): continue    #If already file
            

            f = open(fn,'w')
            ts=0 
            for tag in transcripts.keys():
                if(gene == geneid[tag]):
                    if(ts==maxTran): break #If already recorded maxTran transcripts in gene.fasta file
                    f.write('>' + tag +  '\n')
                    f.write(transcripts[tag]+'\n')
                    ts += 1
            f.close()

        #Now submit Build Super Transcript for each gene in parallel
        print("Now building superTranscript for each multi-transcript gene...")
        jobs = []

        fnames = []
        job_id=1
        for gene in gene_list:
            fname = outdir + '/' + gene + '.fasta'
            fnames.append([fname,"%s of %s" %(job_id,len(gene_list))])
            job_id+=1;

        # BY POOL
        #ncore = 4
#        print("%s clusters will be processed " %len(gene_list) )

        pool = Pool(processes=ncore)
        result = pool.map_async(worker,fnames,chunksize=1)
        pool.close()
        pool.join()
        results = result.get()
        
        #Write Overall Super Duper Tran
        superf = open(outdir + '/' +'SuperDuper.fasta','w')
        supgff = open(outdir + '/' +'SuperDuper.gff','w')

        #Add single cluster genes to the list of results        
        for tag in cluster.keys():
            if(cluster[tag] in single_cluster):
                fnames.append([cluster[tag],''])
                anno=get_annotation_line(cluster[tag],'1',str(len(transcripts[tag])),tag)
                results.append([transcripts[tag],anno,0,1])

        for i,clust in enumerate(fnames):
            #Just use the name of gene, without the preface
            fn = clust[0].split("/")[-1]
            fn = fn.split('.fasta')[0]
            superf.write('>' + fn  + ' NoTrans:' + str(results[i][3]) + ',Whirls:' + str(results[i][2])  + '\n')
            superf.write(results[i][0] + '\n')

        #Write Super gff
        for res in results:
                        supgff.write(res[1])

        Clean(gene_list,outdir,full_clean)

        print("BUILT SUPERTRANSCRIPTS ---- %s seconds ----" %(time.time()-start_time))


def main(args=None):
    #Print Lace Version
    print(" __      __    ____  ____ ")
    print("(  )    / _\  /    )(  __)")
    print("/  (_/\/    \(  (__  ) _) ")
    print("\_____/\_/\_/\_____)(____)")
    print("Lace Version: 1.14")
    

    #Make argument parser
    parser = argparse.ArgumentParser()

    #Add Arguments
    parser.add_argument("TranscriptsFile",help="The name of the fasta file containing all transcripts")
    parser.add_argument("ClusterFile",help="The name of the text file with the transcript to cluster mapping")
    parser.add_argument("--cores",help="The number of cores you wish to run the job on (default = 1)",default=1,type=int)
    parser.add_argument("--alternate","-a",help="Create alternate annotations and create metrics on success of SuperTranscript Building",action='store_true')
    parser.add_argument("--tidy","-t",help="Remove intermediate fasta files after running",action='store_true')
    parser.add_argument("--maxTran",help="Set a maximum for the number of transcripts from a cluster to be included for building the SuperTranscript (default=50).",default=50,type=int)
    parser.add_argument("--outputDir","-o",help="Output Directory",default=".")

    args= parser.parse_args()

    #Make output directory if it currently doesnt exist
    if(args.outputDir):
        if(os.path.exists(args.outputDir)):
            print("Output directory exists")
        else:
            print("Creating output directory")
            os.mkdir(args.outputDir)
    Split(args.TranscriptsFile,args.ClusterFile,args.cores,args.maxTran,args.outputDir,args.tidy)

    if(args.alternate):
        cwd = os.getcwd()

        #Change to output directory
        os.chdir(args.outputDir)
        print("Making Alternate Annotation and checks")
        Checker('SuperDuper.fasta','SuperDuper.gff',args.cores,'SuperFiles')

        #Change back
        os.chdir(cwd)

    print('Done')


if __name__ == '__main__':
    main()
