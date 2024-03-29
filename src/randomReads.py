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
    parser.add_argument("--PFMs",help="Full path(s) to the PFM-file(s). If given, sequences from this PFM are inserted in the middle of the reads (unless another position is specified with --start option). If PFM is not given, random sequences with the given alphabet frequencies are produced.",type=str,default=None,nargs='+')
    parser.add_argument("--seqs",help="Full path to the fasta/fastq-file containing the reads where the PFM is embedded to. If not given, background is generated using the given background alphabet frequencies.",type=str,default=None)
    parser.add_argument("--L",help="Length of the output sequences (default=100, overruled by the sequence length of --seqs if given).",type=int,default=100)
    parser.add_argument("--N",help="Number of sequences generated (default=100000, overruled by the number of sequences in --seqs if given).",type=int,default=100000)
    parser.add_argument("--distance",help="If given, a pair of PFMs is always inserted at a fixed istance from each other (default=None)",type=int,default=None)
    parser.add_argument("--starts",help="Start position(s) of the inserted PFMs (default=random). If multiple PFMs given, each needs to be given its own start position.",type=int,default=None,nargs='+')
    parser.add_argument("--bgfreqs",help="Background alphabet frequencies, default is flat background distribution. Order for DNA: A, C, G, T. Order for protein: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y. Order for RNA: A, C, G, U.",type=float,default=None,nargs='+')
    parser.add_argument("--concensus",help="If yes, always insert the concensus of the PWM(s). If no (=defaults), sample from the PFM(s).",type=str,choices=['yes','no'],default='no')
    parser.add_argument("--addToReadName",help="String added to read names to distinguish them from background reads (default=embed).",type=str,default=":embed")
    parser.add_argument("--alphabet",help="Alphabet used, choices are DNA (=default) or protein.",type=str,choices=['DNA','protein','RNA'],default='DNA')
    parser.add_argument("--seed",help="Seed for the random number generator (default=42).",type=int,default=42)
    
    args = parser.parse_args()
    np.random.seed(args.seed)
    
    #reading in the pfm/pwm if it is given
    if args.PFMs!=None:
        PFMs = [] #list containing the normalised PFM matrices
        Ls = [] #list containing the PFM lengths
        for PFMfile in args.PFMs:
            with open(PFMfile,'rt') as csvfile:
                r = csv.reader(csvfile,delimiter='\t')
                PFMs.append([])
                for row in r: PFMs[-1].append([float(i) for i in row])
            PFMs[-1] = np.array(PFMs[-1])
            Ls.append(PFMs[-1].shape[1])
            if args.concensus=='no':
                #normalizing the matrix entries so that each column sums up to one
                for i in range(0,Ls[-1]): PFMs[-1][:,i] /= sum(PFMs[-1][:,i])
            else:
                #setting all other matrix entries but the concensus sequence as 0
                concensus = np.argmax(PFMs[-1],axis=0)
                PFMs[-1] = np.zeros(shape=PFMs[-1].shape)
                for i in range(0,PFMs[-1].shape[1]): PFMs[-1][concensus[i],i] = 1.0

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
                    if args.starts==None:
                        if len(PFMs)>1:
                            #print("Define start positions for the embedded matches!")
                            #exit
                            #we simply draw start positions by random until they don't overlap
                            #this works currently up to two motifs only
                            while True:
                                starts = [np.random.randint(0,high=len(seq)-Ls[0]),np.random.randint(0,high=len(seq)-Ls[1])]
                                if starts[0]<starts[1]:
                                    if starts[0]+Ls[0]<starts[1]: break
                                else:
                                    if starts[1]+Ls[1]<starts[0]: break
                                    
                        else: starts = [np.random.randint(0,high=len(seq)-Ls[0])]
                    else: starts = args.starts
                    first = False

                if args.starts==None:
                    if len(PFMs)<2: starts = [np.random.randint(0,high=len(seq)-Ls[0])]
                    else:
                        if args.distance!=None:
                            #this only works for pairs currently
                            #first draw position for the first of pair
                            first_tf_ind = np.random.randint(0,high=2)
                            starts = [0,0]
                            #print(first_tf_ind)
                            if first_tf_ind==1:
                                starts[1] = np.random.randint(0,high=len(seq)-Ls[0]-Ls[1]-args.distance)
                                starts[0] = starts[1]+Ls[1]+args.distance
                            else:
                                starts[0] = np.random.randint(0,high=len(seq)-Ls[0]-Ls[1]-args.distance)
                                starts[1] = starts[0]+Ls[0]+args.distance
                        else:
                            while True:
                                starts = [np.random.randint(0,high=len(seq)-Ls[0]),np.random.randint(0,high=len(seq)-Ls[1])]
                                if starts[0]<starts[1]:
                                    if starts[0]+Ls[0]<starts[1]: break
                                else:
                                    if starts[1]+Ls[1]<starts[0]: break

                header = fasta.id
                header = ">"+header+args.addToReadName
                newseq = seq
                #print("starts:"+str(starts))
                for p in range(0,len(PFMs)):
                    randoms = np.random.rand(Ls[p])
                    PFM_seq = ""
                    for i in range(0,len(randoms)):
                        nucl_index = 0
                        count = PFMs[p][0,i]
                        while True:
                            if count>=randoms[i]:
                                if nucl_index==0: PFM_seq += 'a'
                                elif nucl_index==1: PFM_seq += 'c'
                                elif nucl_index==2: PFM_seq += 'g'
                                elif nucl_index==3: PFM_seq += 't'
                                break
                            nucl_index += 1
                            count += PFMs[p][nucl_index,i]
                    newseq = newseq[:starts[p]]+PFM_seq+newseq[starts[p]+Ls[p]:]
                #saving the new sequence
                w.writerow([header])
                w.writerow([newseq])
    else:
        #generating the background sequence from the the given alphabet frequencies
        #creating one args.L length PWM where all sequences are drawn from
        if args.alphabet=='DNA' or args.alphabet=='RNA':
            A = 4 #alphabet length
        elif args.alphabet=='protein':
            A = 20 #alphabet length
        full_PFM = np.ones(shape=(A,args.L))
        if args.bgfreqs!=None:
            if len(args.bgfreqs)!=A:
                print("Wrong number of background frequencies given!")
                exit
            for i in range(0,A): full_PFM[i,:] *= args.bgfreqs[i]
        else:
            for i in range(0,A): full_PFM[i,:] *= (1.0/A)
            
        #inserting the PFM-counts if needed
        if args.PFMs!=None:
            for p in range(0,len(args.PFMs)):
                #determining the insertion positions for the PFM-derived sequences
                if p==0:
                    if args.starts==None:
                        if len(PFMs)>1:
                            print("Define start positions for the embedded matches!")
                            exit
                        starts = [int(np.random.randint(0,high=len(seq)-Ls[p]))]
                    else: starts = args.starts

                #adding the PFM entries
                full_PFM[:,starts[p]:starts[p]+Ls[p]] = PFMs[p]


                
        #creating the sequences
        with open(args.outfile,'wt') as outfile:
            w = csv.writer(outfile,delimiter='\t')
            for i in range(0,args.N):
                header = ">"+str(i)
                randoms = np.random.rand(args.L)
                PFM_seq = ""
                #print(randoms)
                for i in range(0,len(randoms)):
                    index = 0
                    count = full_PFM[0,i]
                    while True:
                        #print("count="+str(count))
                        #print("random="+str(randoms[i]))
                        if count>=randoms[i]:
                            if args.alphabet=='DNA':
                                if index==0: PFM_seq += 'A'
                                elif index==1: PFM_seq += 'C'
                                elif index==2: PFM_seq += 'G'
                                elif index==3: PFM_seq += 'T'
                                break
                            elif args.alphabet=='RNA':
                                if index==0: PFM_seq += 'A'
                                elif index==1: PFM_seq += 'C'
                                elif index==2: PFM_seq += 'G'
                                elif index==3: PFM_seq += 'U'
                                break
                            elif args.alphabet=='protein':
                                if index==0: PFM_seq += 'A'
                                elif index==1: PFM_seq += 'C'
                                elif index==2: PFM_seq += 'D'
                                elif index==3: PFM_seq += 'E'
                                elif index==4: PFM_seq += 'F'
                                elif index==5: PFM_seq += 'G'
                                elif index==6: PFM_seq += 'H'
                                elif index==7: PFM_seq += 'I'
                                elif index==8: PFM_seq += 'K'
                                elif index==9: PFM_seq += 'L'
                                elif index==10: PFM_seq += 'M'
                                elif index==11: PFM_seq += 'N'
                                elif index==12: PFM_seq += 'P'
                                elif index==13: PFM_seq += 'Q'
                                elif index==14: PFM_seq += 'R'
                                elif index==15: PFM_seq += 'S'
                                elif index==16: PFM_seq += 'T'
                                elif index==17: PFM_seq += 'V'
                                elif index==18: PFM_seq += 'W'
                                elif index==19: PFM_seq += 'Y'
                                break
                        index += 1
                        count += full_PFM[index,i]
                    #print(i)
                    #print(PFM_seq)
                #saving the new sequence
                w.writerow([header])
                w.writerow([PFM_seq])
            
#end

randomReads()
