#! /usr/bin/python
"""
Created on Tue Nov 24 12:22:10 2015

@author: charles
The purpose of this program is to locate putative transcription start sites from an alignment of 5'-capture reads against
a reference in SAM format. It counts the first base of the first read of the paired end reads to be the TSS
"""
import sys
import numpy as np
import collections
import matplotlib.pyplot as plt
def main(args):
    SAM_File = open(sys.argv[1])
    CHR_array_forward = {}
    CHR_array_reverse = {}
    for line in SAM_File:
        if line[0:3] == '@SQ':
            line_array = line.strip().split("\t")
            header = line_array[1][3::]
            CHR_array_forward[header] = collections.defaultdict(int)
            CHR_array_reverse[header] = collections.defaultdict(int)
        elif line[0] != '@':
            line_array=line.split("\t")
            bit_flag = int(line_array[1])
            CIGAR = str(line_array[5])
            chromosome = str(line_array[2])
            left_most_base = int(line_array[3])
            if (bit_flag & 64 or not bit_flag & 1) and not (bit_flag & 4):
                #print("First read: ",bit_flag,distance)
                displacement = CIGAR_length(CIGAR)
                if bit_flag & 16:
                    CHR_array_reverse[chromosome][left_most_base + displacement]+= 1
                else:
                    CHR_array_forward[chromosome][left_most_base]+= 1
    print("track type=bedGraph")
    for peak in find_peaks(CHR_array_forward):
        print(peak[0]+"\t"+str(peak[1] - 1)+"\t"+str(peak[1])+"\t"+str(peak[2]))
    for peak in find_peaks(CHR_array_reverse):
        print(peak[0]+"\t"+str(peak[1] - 1)+"\t"+str(peak[1] - 2)+"\t"+str(-peak[2]))    


def find_peaks(CHR_array):
    for key in CHR_array.keys():
        chromosome = CHR_array[key]
        for i in chromosome.keys():
            if chromosome[i] >= 5 and chromosome[i] >= chromosome[i - 1] and chromosome[i] >= chromosome[i + 1]:
                yield(key,i,chromosome[i])

def CIGAR_length(S):
    length = 0
    current_bases = ""
    for C in S:
        if C in "ATGC=X":
            length += 1
        elif C in "MDN":
            length += int(current_bases)
            current_bases = ""
        elif C in "SHI":
            current_bases = ""
        elif C in "0123456789":
            current_bases+=(C)
    return(length)
if (__name__ == "__main__"):
    sys.exit(main(sys.argv))