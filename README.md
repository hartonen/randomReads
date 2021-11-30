
# randomReads

## INTRODUCTION

RandomReads is a simple python script for embedding transcription factor binding motifs (input as PFM matrices) to given DNA sequences or to a given mononucleotide background. Python packages needed (and the versions tested with):

```
biopython==1.72
numpy==1.15.1
```


No installation is required. To use the randomReads.py script the src/ directory should be added to $PATH or the .py-file from the src/ directory should be copied to a directory already listed in your $PATH. If we assume you have unpacked randomReads to directory home/randomReads/, you should add the following line to your ~/.bashrc or whatever file you use to set your environment variables

`export PATH=$PATH:home/randomReads/src`

In examples below we use [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html) to fetch random background sequences from the human genome.

## USAGE

```
usage: randomReads.py [-h] [--PFMs PFMS [PFMS ...]] [--seqs SEQS] [--L L]
                      [--N N] [--distance DISTANCE]
                      [--starts STARTS [STARTS ...]]
                      [--bgfreqs BGFREQS [BGFREQS ...]] [--concensus {yes,no}]
                      [--addToReadName ADDTOREADNAME]
                      [--alphabet {DNA,protein,RNA}] [--seed SEED]
                      outfile

positional arguments:
  outfile               Output fasta-file name.

optional arguments:
  -h, --help            show this help message and exit
  --PFMs PFMS [PFMS ...]
                        Full path(s) to the PFM-file(s). If given, sequences
                        from this PFM are inserted in the middle of the reads
                        (unless another position is specified with --start
                        option). If PFM is not given, random sequences with
                        the given alphabet frequencies are produced.
  --seqs SEQS           Full path to the fasta/fastq-file containing the reads
                        where the PFM is embedded to. If not given, background
                        is generated using the given background alphabet
                        frequencies.
  --L L                 Length of the output sequences (default=100, overruled
                        by the sequence length of --seqs if given).
  --N N                 Number of sequences generated (default=100000,
                        overruled by the number of sequences in --seqs if
                        given).
  --distance DISTANCE   If given, a pair of PFMs is always inserted at a fixed
                        istance from each other (default=None)
  --starts STARTS [STARTS ...]
                        Start position(s) of the inserted PFMs
                        (default=random). If multiple PFMs given, each needs
                        to be given its own start position.
  --bgfreqs BGFREQS [BGFREQS ...]
                        Background alphabet frequencies, default is flat
                        background distribution. Order for DNA: A, C, G, T.
                        Order for protein: A, C, D, E, F, G, H, I, K, L, M, N,
                        P, Q, R, S, T, V, W, Y. Order for RNA: A, C, G, U.
  --concensus {yes,no}  If yes, always insert the concensus of the PWM(s). If
                        no (=defaults), sample from the PFM(s).
  --addToReadName ADDTOREADNAME
                        String added to read names to distinguish them from
                        background reads (default=embed).
  --alphabet {DNA,protein,RNA}
                        Alphabet used, choices are DNA (=default) or protein.
  --seed SEED           Seed for the random number generator (default=42).															      
```


## EXAMPLE 1: EMBEDDING CTCF MOTIF TO GENOMIC BACKGROUND

We use bedtools to fetch random regions from the hg19 human genome assembly. First we extract 1000000 random 100 bp long regions from the genome:

`bedtools random -l 100 -n 1000000 -g path/hg19.chrom.sizes > hg19_bg.bed`

Then we fetch the underlying sequences:

`bedtools getfasta -fi hg19.fasta -bed hg19_bg.bed -fo hg19_bg.fasta`

We use these sequences as background templates and place the CTCF motif in the middle of the 100 bp sequences:

`randomReads.py hg19_bg_CTCF_embedded.fasta --PFM CTCF_AJ_TAGCGA20NGCT_NGCGCCMYCTAGYGGTN_m2_c4_Cell2013.pfm --seqs hg19_bg.fasta`

The sequences with the CTCF motif embedded in the middle are now in file hg19_bg_CTCF_embedded.fasta.
