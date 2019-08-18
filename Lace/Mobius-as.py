#!/usr/bin/env python
#Author: Anthony Hawkins
#A little script to construct an annotation from an SJ.out.tab file, a standard STAR output. 
# This creates the dynamic block annotation

#Edited by Damayanthi to read the transcript start and end information from the original annotation file.
import argparse
import pandas as pd
import os, sys
import bisect

def Mobius(sjfile,gfile,ft,ann_trans,read_thresh,flat_ann,outputfileName):

    #First read in the SJ.out.tab as sys.argv[1]
    sj = pd.read_csv(sjfile,sep='\t',header=None,names=['Gene','Start','End','strand','intron motif','Annotated','Unique','Multi','Overhang'])
    #Read in SuperTranscript fasta file
    sf = open(sys.argv[2],'r')
    glength = {} #A dictionary holding the SuperTranscript gene lengths
    gene=''
    for line in sf:
        if('>' in line):
            gene= (line.split()[0]).split('>')[1]
            glength[gene] = ''
        else:
            glength[gene] = glength[gene] + line.split('\n')[0].split('\r')[0]

    #Create gtf file
    
    gtf = open(outputfileName,'w')
    slist = {}
    #Make a dictionary for each gene, which holds a list of splice junction start and end points
    #For each row
    for i in range(0,len(sj['Gene'])):
        curr_gene = sj.iloc[i,0]
        if(curr_gene not in slist.keys()): slist[curr_gene] = [1]
        if(sj.iloc[i,6] > read_thresh): #More than 5 reads (only uniquely mapping reads  spanning junction)
            slist[curr_gene].append(int(sj.iloc[i,1])) #This is actually the intron start part (i.e one base over)
            slist[curr_gene].append(int(sj.iloc[i,2])+1) #This is the end of the exon-exon junction, so need to add one for the actual exon start place

    #Read transcript annotation 
    #If forcing transcript start and ends
    igtf = pd.read_csv(flat_ann,sep='\t',header=None,names=['Gene','Source','Type','Start','End','v1','v2','v3','v4'])
    #Create a dict in the form {gene_id:transcript_start_positions}
    defaultistTranscriptStarts={} #dict to hold start sites of transcripts of a given gene
    defaultistTranscriptEnds={} #dict to hold end sites of transcripts of a given gene.


    prev_gene= ""
    next_gene=""
    i=0
    prev_trans= igtf.iloc[0,8].split(';')[1].split('"')[1]
    prev_end= 0
    for i in range(0,len(igtf['Gene'])-1):
        curr_gene = igtf.iloc[i,0]
        next_gene = igtf.iloc[i+1,0]
        curr_trans = igtf.iloc[i,8].split(';')[1].split('"')[1]
        #If a new transcript add start site

        if (curr_gene!=prev_gene): #add start new gene
            if curr_gene not in  defaultistTranscriptStarts.keys(): defaultistTranscriptStarts[curr_gene]= [int(igtf.iloc[i,3])]
            else: defaultistTranscriptStarts[curr_gene].append(int(igtf.iloc[i,3])) #Transcript start position.

        if (curr_gene!=next_gene): #Add transcript end sites.
            if curr_gene not in defaultistTranscriptEnds.keys():defaultistTranscriptEnds[curr_gene] = [int(igtf.iloc[i,4])]
            else: defaultistTranscriptEnds[curr_gene].append(int(igtf.iloc[i,4]))


        if (((curr_gene==prev_gene) | (curr_gene ==next_gene)) & (curr_trans!=prev_trans)): #Change of transcripts of a same gene
            defaultistTranscriptStarts[curr_gene].append(int(igtf.iloc[i,3])) #Transcript start position.

            if (curr_gene ==prev_gene): # Add end sites
                if curr_gene not in defaultistTranscriptEnds.keys():defaultistTranscriptEnds[curr_gene] = [int(prev_end)]
                else: defaultistTranscriptEnds[curr_gene].append(int(prev_end))

        prev_gene = curr_gene
        prev_end = igtf.iloc[i,4]
        prev_trans= curr_trans

    #Account for last entry
    curr_gene = igtf.iloc[i,0]

    if (curr_gene!=next_gene):
        if curr_gene not in  defaultistTranscriptStarts.keys(): defaultistTranscriptStarts[curr_gene]= [int(igtf.iloc[i+1,3])]
        else: defaultistTranscriptStarts[curr_gene].append(int(igtf.iloc[i,3])) #Transcript start position.
        if curr_gene not in defaultistTranscriptEnds.keys():defaultistTranscriptEnds[curr_gene] = [int(igtf.iloc[i+1,4])]
        else: defaultistTranscriptEnds[curr_gene].append(int(igtf.iloc[i,4]))

    #Keep unique elements and sort the transcript start/end sites dict.

    for key in defaultistTranscriptStarts.keys():
        defaultistTranscriptStarts[key] = list(set(defaultistTranscriptStarts[key]))
        defaultistTranscriptStarts[key].sort()

        defaultistTranscriptEnds[key] = list(set(defaultistTranscriptEnds[key]))
        defaultistTranscriptEnds[key].sort()



    for key in glength.keys():
        exon_counter = 1
        if(key in slist.keys()):
            slist[key] = list(set(slist[key]))
            slist[key].sort()

            #Add blocks to account for alternative start and end events
            #Retrieve first exon


            if(key in defaultistTranscriptStarts.keys()):

                #add all the transcript start sites before first exon in splice junction list

                if len(slist[key])>1: #Avoid cases where there were no significant splice junctions
                    starts =  [i for i in defaultistTranscriptStarts[key] if i <= slist[key][1]]
                    for j in starts: slist[key].append(j)

                    ends = [i for i in defaultistTranscriptEnds[key] if i >= slist[key][-1]]
                    for j in ends: slist[key].append(j)


            #Remove if lenght of ST is already added to slist, to remove confusion when writing the last entry corresponding to a gene.
            slist[key] = [i for i in slist[key] if i != len(glength[key])]

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


def main(args=None):
	#ASCII art
        print(" ___ ___   ___   ____   ____  __ __  _____")
        print("|   |   | /   \ |    \ |    ||  |  |/  __/")
        print("| _   _ ||     ||  o  ) |  | |  |  (   \_ ")
        print("|  \_/  ||  o  ||     | |  | |  |  |\__  |")
        print("|   |   ||     ||  O  | |  | |  :  |/  \ |")
        print("|   |   ||     ||     | |  | |     |\    |")
        print("|_______| \___/ |_____||____| \__,_| \___|")


        #Make argument parser
        parser = argparse.ArgumentParser()
        #Add Arguments
        parser.add_argument("SpliceJunctions",help="The name of the Splice Junctions tab file (in the format of the one STAR produces)")
        parser.add_argument("GenomeFasta",help="A fasta file containing the sequence for all genes in genome")
        parser.add_argument("-forceTrans","-ft",help="Force blocks where annotated transcripts start and end",action='store_true')
        parser.add_argument("-AnnoTrans","-a",help="Flattened SuperTranscript annotation file",default="")
        parser.add_argument("-readThresh","-reads",help="The minimum number of reads in all combined samples required to support a splice junction (default=5)",default=5)
	parser.add_argument("flat_ann",help="Annotation to consider")
	parser.add_argument("outputfileName",help="output file name")
        args= parser.parse_args()

        print('Constructing Dynamic Blocks- AS')
        Mobius(args.SpliceJunctions,args.GenomeFasta,args.forceTrans,args.AnnoTrans,args.readThresh,args.flat_ann,args.outputfileName)
        print('Done')


if __name__ == '__main__':
	main()
