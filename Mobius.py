#Author: Anthony Hawkins
#A little script to construct an annotation from an SJ.out.tab file, a standard STAR output. 
# This creates the dynamic block annotation

import argparse
import pandas as pd
import os, sys

def Mobius(sjfile,gfile,ft,flat_ann):

	#First read in the SJ.out.tab as sys.argv[1]
    sj = pd.read_csv(sjfile,sep='\t',header=None,names=['Gene','Start','End','strand','intron motif','Annotated','Unique','Multi','Overhang'])

    #Read in SuperTranscript fasta file
    sf = open(sys.argv[2],'r')
    glength = {} #A dictionary holding the SuperTranscript gene lengths
    gene=''
    for line in sf:
        if('>' in line):
            gene= (line.split(' ')[0]).split('>')[1]
            glength[gene] = ''
        else:
            glength[gene] = glength[gene] + line.split('\n')[0].split('\r')[0]

    #Create gtf file
    gtf = open('Spliced.gtf','w')
    slist = {}


    #Make a dictionary for each gene, which holds a list of splice junction start and end points
    #For each row
    for i in range(0,len(sj['Gene'])):
        curr_gene = sj.iloc[i,0]
        if(curr_gene not in slist.keys()): slist[curr_gene] = [1]
        if((sj.iloc[i,7] + sj.iloc[i,8]) > 5): #More than 5 reads (either unique or multi spanning junction)
            slist[curr_gene].append(int(sj.iloc[i,1])) #This is actually the intron start part (i.e one base over)
            slist[curr_gene].append(int(sj.iloc[i,2])+1) #This is the end of the exon-exon junction, so need to add one for the actual exon start place

    #If forcing transcript start and ends
    if(ft): 
        #Make a dictionary with the transcript list per gene
        igtf = pd.read_csv(flat_ann,sep='\t',header=None,names=['Gene','Source','Type','Start','End','v1','v2','v3','v4'],skiprows=2)
        prev_gene =''
        prev_trans=''
        for i in range(0,len(igtf)):
            curr_gene = igtf.iloc[i,0]
            curr_trans = igtf.iloc[i,8].split(';')[1].split('"')[1]

	    #If gene not already in spliced dict, then make a key with an empty list value
            if (curr_gene not in slist.keys()): slist[curr_gene]=[]

            #Switch transcripts in gene
            if(curr_trans != prev_trans and curr_gene==prev_gene): #Just switched from one transcript in gene to another one in same gene
                slist[curr_gene].append(prev_end) #End of previous transcript
                slist[curr_gene].append(igtf.iloc[i,3]) #Start of new transcript

            #Switch gene
            if(curr_gene != prev_gene): #Beginning a whole new gene
                slist[curr_gene].append(igtf.iloc[i,3]) #Start site of first transcript in new gene
                if(i != 0): slist[prev_gene].append(prev_end) #Make sure add in the end point of the last transcript on the last gene
		
            prev_gene = curr_gene
            prev_trans = curr_trans
            prev_end = igtf.iloc[i,4]

     #Now sort each list for each gene, only keep unique elements
    for key in glength.keys():
        exon_counter = 1
        if(key in slist.keys()):
            slist[key] = list(set(slist[key]))
            slist[key].sort()

	    #Now for each coord in list make a block
            for i in range(1,len(slist[key])):
                    ann = str(key) + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(slist[key][i-1])  + '\t' + str(slist[key][i]-1) + '\t' + '.' + '\t' + '.' + '\t' + '0' + '\t' + 'gene_id "' +str(key) +'"; transcript_id "' + str(key) + '"; exon_number ' + str(exon_counter) + '; exon_id "' + str(key)+':'+str(exon_counter)+ '"' + '\n' #Note: need to minus one off end to account for the fact that the exon ends before the exon-exon boundary exists
                    exon_counter +=1
                    gtf.write(ann)


        #For the list splice junnction, we need to make a block from the last sj to the end of the ST
        if(key not in slist): last = 1
        else: last = slist[key][-1]
 
	#If not forcing transcript start and ends
        if(not ft):
            if(last != len(glength[key])): 
                ann = str(key) + '\t' + 'SuperTranscript' + '\t' + 'exon' + '\t' + str(last)  + '\t' + str(len(glength[key])) + '\t' + '.' + '\t' + '.' + '\t' + '0' + '\t' + 'gene_id "' +str(key) +'"; transcript_id "' + str(key) + '"; exon_number ' + str(exon_counter) + '; exon_id "' + str(key)+':'+str(exon_counter)+ '"' + '\n'
                gtf.write(ann)


if __name__ == '__main__':

        #Make argument parser
        parser = argparse.ArgumentParser()

        #Add Arguments
        parser.add_argument("SpliceJunctions",help="The name of the Splice Junctions tab file (in the format of the one STAR produces)")
        parser.add_argument("GenomeFasta",help="A fasta file containing the sequence for all genes in genome")
        parser.add_argument("-forceTrans","-ft",help="Force blocks where annotated transcripts start and end",action='store_true')
        parser.add_argument("-AnnoTrans","-a",help="Flattened SuperTranscript annotation file",default="")
        args= parser.parse_args()

        print('Constructing Dynamic Blocks')
        Mobius(args.SpliceJunctions,args.GenomeFasta,args.forceTrans,args.AnnoTrans)
        print('Done')
