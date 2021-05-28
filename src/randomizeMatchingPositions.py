#!/usr/bin/env python

import argparse
import pyfastx
import csv

import numpy as np

def randomizeMatchingPositions():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #MANDATORY PARAMETERS
    parser.add_argument("outfile",help="Output fasta-file name.",type=str)

    #OPTIONAL PARAMETERS
    parser.add_argument("--wt",help="Full path to a fasta-file containing the wild type sequence.",type=str)
    parser.add_argument("--seqs",help="Full path to the fasta-file containing the sequences where we want to randomize the positions matching to wild type.",type=str)
    parser.add_argument("--N",help="Exact number of mismatches in seqs needed for including to output.",type=int,default=2)
    parser.add_argument("--addToReadName",help="String added to read names to distinguish them from input reads (default=:randomized).",type=str,default=":randomized")
    parser.add_argument("--alphabet",help="Alphabet used as a string containing each possible character (case sensitive, default=ACGT).",type=str,default='ACGT')

    args = parser.parse_args()

    #read in the wild type sequence
    for name,seq in pyfastx.Fasta(args.wt):
        wtseq = seq

    #read in rest of the sequences, save to outfile those that have N mismatches to wtseq
    #and randomize other positions from them
    with open(args.outfile,'wt') as outfile:
        w = csv.writer(outfile,delimiter='\t')
        for name,seq in pyfastx.Fasta(args.seqs):
            rands = np.random.randint(0,high=len(args.alphabet),size=len(seq)) #draw the random sequence
            newseq = ""
            N_mismatch = 0 #mismatch counter
            for i in range(0,len(seq)):
                if seq[i]!=wtseq[i]:
                    N_mismatch += 1
                    newseq += seq[i]
                else: newseq += args.alphabet[rands[i]]
                if N_mismatch>args.N: break
            if n_mismatch==args.N:
                #save the sequence
                w.writerow(['>'+name+args.addToReadName])
                w.writerow([newseq])
#end

randomizeMatchingPositions()
