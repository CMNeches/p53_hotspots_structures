#!/usr/bin/env python3

#import packages
from Bio import SeqIO
from Bio.Seq import Seq
import os.path
import string
import argparse

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("isoform_name", help = "isoform name")
args = parser.parse_args()
isoform_name = args.isoform_name

#read in wt nucleotide sequence
wt_seq = list(SeqIO.parse("../../data/sequences/TP53_coding_sequence.fasta", "fasta"))

#isoform_names = ["isoform_a", "isoform_b", "isoform_c", "isoform_d", "isoform_e", "isoform_f", "isoform_g", "isoform_h", "isoform_i", "isoform_j", "isoform_k", "isoform_l", "isoform_m", "isoform_n", "isoform_o"]

isoform_names = ["isoform_" + l for l in list(string.ascii_lowercase)[:15]]

wt_seq = wt_seq[isoform_names.index(isoform_name)].seq
aa_wt_seq = wt_seq.translate()

ref_nucs = ["A", "C", "G", "T"]

for i, n1 in enumerate(wt_seq): #enumerate starts counting at zero
    for  n2 in ref_nucs:
        mut_seq = wt_seq[:i] + n2 + wt_seq[i+1:] #+1 leaves a hole for n2
        mut_seq = str(mut_seq.translate()) #it is going from 3' to 5'
       
        #looking for amino acid change
        count = 0
        for a1, a2 in zip(aa_wt_seq, mut_seq):
            count +=1
            if a1 != a2:
                pos = count
                break

        #don't do anything if the sequences are the same
        if aa_wt_seq == mut_seq:
            pass
        else:
            #stop at the stop codon
            if "*" in mut_seq:
                end = mut_seq.index("*")
                mut_seq = mut_seq[:end]
            else:
                pass
            
            #write to a fasta file
            mut = "{}{}{}".format(a1, pos, a2).replace("*", "STOP")
            name = ">p.{}".format(mut)
            file_path = "../../results/mutant_sequences/SNV/{}/{}.fasta".format(isoform_name, mut)
            if os.path.exists(file_path) == False:
                with open(file_path, "w") as f:
                    f.write(name + "\n" + str(mut_seq) + "\n")
            else:
                pass
