#!/usr/bin/env python

import argparse
from Bio import SeqIO
import csv
import gzip
import numpy as np

def randomReads():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #MANDATORY PARAMETERS
    parser.add_argument("outfile",help="Output fasta-file name.",type=str)

    #OPTIONAL PARAMETERS
    parser.add_argument("--PFM",help="Full path to the PFM-file. If given, sequences from this PFM are inserted in the middle of the reads (unless another position is specified with --start option). If PFM is not given, random sequences with the given mononucleotide frequencies are produced.",type=str,default=None)
    parser.add_argument("--seqs",help="Full path to the fasta/fastq-file containing the reads where the PFM is embedded to. If not given, background is generated using the given mononucleotide frequencies.",type=str,default=None)
    parser.add_argument("--L",help="Length of the output sequences (default=100, overruled by the sequence length of --seqs if given).",type=int,default=100)
    parser.add_argument("--N",help="Number of sequences generated (default=100000, overruled by the number of sequences in --seqs if given).",type=int,default=100000)
    parser.add_argument("--start",help="Start position of the inserted PFM-sequences (default=middle of the read).",type=int,default=None)
    parser.add_argument("--bgfreqs",help="Background nucleotide frequencies: A, C, G, T (default=0.25,0.25,0.25,0.25).",type=float,default=[0.25,0.25,0.25,0.25],nargs=4)
    
    args = parser.parse_args()

    #reading in the pfm/pwm if it is given
    if args.PFM!=None:
        with open(args.PFM,'rt') as csvfile:
            r = csv.reader(csvfile,delimiter='\t')
            PFM = []
            for row in r: PFM.append([float(i) for i in row])
        PFM = np.array(PFM)
        #normalizing the matrix entries so that each column sums up to one
        L = PFM.shape[1]
        for i in range(0,L): PFM[:,i] /= sum(PFM[:,i])
        #print(PFM)
        #print(PFM.shape)

    #either reading in or generating the output sequences one by one
    if args.seqs!=None:
        #reading in the background sequences from file
        ftype = args.seqs.split('.')[-1]
        fasta_sequences = SeqIO.parse(open(args.seqs),ftype)

        with open(args.outfile,'wt') as outfile:
            w = csv.writer(outfile,delimiter='\t')
            first = True
            
            for fasta in fasta_sequences:
                seq = str(fasta.seq).upper()
                if first:
                    #determining the insertion position for the PFM-derived sequences
                    if args.start==None: start = int(len(seq)/2-L/2)
                    else: start = args.start
                    first = False
                    
                header = ">"+str(i)
                randoms = np.random.rand(L)
                PFM_seq = ""
                for i in range(0,len(randoms)):
                    nucl_index = 0
                    count = PFM[0,i]
                    while True:
                        if count>=randoms[i]:
                            if nucl_index==0: PFM_seq += 'A'
                            elif nucl_index==1: PFM_seq += 'C'
                            elif nucl_index==2: PFM_seq += 'G'
                            elif nucl_index==3: PFM_seq += 'T'
                            break
                        nucl_index += 1
                        count += PFM[nucl_index,i]
                #saving the new sequence
                w.writerow([header])
                w.writerow([seq[:start]+PFM_seq+seq[start+L:]])
    else:
        #generating the background sequence from the the given nucleotide frequencies
        #creating one args.L length PWM where all sequences are drawn from
        full_PFM = np.ones(shape=(4,args.L))
        for i in range(0,4): full_PFM[i,:] *= args.bgfreqs[i]
        #inserting the PFM-counts if needed
        if args.PFM!=None:
            #determining the insertion position for the PFM-derived sequences
            if args.start==None: start = int(args.L/2-L/2)
            else: start = args.start

            #adding the PFM entries
            full_PFM[:,start:start+L] = PFM


        print(full_PFM)
        print(full_PFM.shape)
        
        #creating the sequences
        with open(args.outfile,'wt') as outfile:
            w = csv.writer(outfile,delimiter='\t')
            for i in range(0,args.N):
                header = ">"+str(i)
                randoms = np.random.rand(args.L)
                PFM_seq = ""
                #print(randoms)
                for i in range(0,len(randoms)):
                    nucl_index = 0
                    count = full_PFM[0,i]
                    while True:
                        #print("count="+str(count))
                        #print("random="+str(randoms[i]))
                        if count>=randoms[i]:
                            if nucl_index==0: PFM_seq += 'A'
                            elif nucl_index==1: PFM_seq += 'C'
                            elif nucl_index==2: PFM_seq += 'G'
                            elif nucl_index==3: PFM_seq += 'T'
                            break
                        nucl_index += 1
                        count += full_PFM[nucl_index,i]
                    #print(i)
                    #print(PFM_seq)
                #saving the new sequence
                w.writerow([header])
                w.writerow([PFM_seq])
            
#end

randomReads()
