#! /usr/bin/python
"""
Created on Mon Jan 11 11:16:34 2016

@author: charles
The purpose of this program is to calculate the minimum distance between features in two bed files.
This program was created so that I can get a sense for how closely two sets of peaks overlap.
This program uses an BST data structure to store and search points
"""
import sys
import collections
from scipy import spatial
from matplotlib import pyplot
import math
def main(args):
    Reference_hash = collections.defaultdict(default_ref_hash)
    Tree_hash = collections.defaultdict(default_ref_hash)
    Reference_file = open(sys.argv[1])
    Reference_file.readline()
    Distance_Array = []
    for line in Reference_file:
        line_array = line.split("\t")
        sign = "+"
        if line_array[3][0] == "-":
            sign = "-"        
        Reference_hash[line_array[0]][sign].append([int(line_array[1])])
        
    for i in Reference_hash.keys():
        for j in Reference_hash[i].keys():
            if len(Reference_hash[i][j]) > 0:
                Tree_hash[i][j] = spatial.KDTree(Reference_hash[i][j])
    Query_file = open(sys.argv[2])
    Query_file.readline()
    for line in Query_file:
        line_array = line.split("\t")
        sign = "+"
        if line_array[3][0] == "-":
            sign = "-"
        if Tree_hash[line_array[0]][sign] != []:
            closest_point = Tree_hash[line_array[0]][sign].query([[int(line_array[1])]])
            #print(Reference_hash[line_array[0]][sign][closest_point[1][0]][0],int(line_array[1]))
            if Reference_hash[line_array[0]][sign][closest_point[1][0]][0] > int(line_array[1]) and sign == "+":
                print(closest_point[0][0]*-1)
            elif Reference_hash[line_array[0]][sign][closest_point[1][0]][0] < int(line_array[1]) and sign == "-":
                print(closest_point[0][0]*-1)
            else:
                print(closest_point[0][0])
    
def default_ref_hash():
    temp_hash = {}
    temp_hash["-"] = []
    temp_hash['+'] = []
    return temp_hash
if (__name__ == "__main__"):
    sys.exit(main(sys.argv))