from __future__ import division
__author__ = 'alipirani'
import statistics
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np

def parser():
    parser = argparse.ArgumentParser(description='Estimate ptr from coverage bed file for qPCR loci')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-bedfile', action='store', dest="bedfile", help='Coverage Bed file', required=True)
    required.add_argument('-outfile', action='store', dest="outfile", help='PTR for qPCR loci', required=True)
    required.add_argument('-OriC_coordinates', action='store', dest="OriC_coordinates", help='OriC Start and End coordinates', required=True)
    required.add_argument('-ter_coordinates', action='store', dest="ter_coordinates", help='ter Start and End coordinates', required=True)
    return parser

# Main Method
# The algorithm follows the procedure as described in this publication: 
# Growth dynamics of gut microbiota in health and disease inferred from single metagenomic samples 
# http://science.sciencemag.org/content/349/6252/1101.long with a few minor changes.
# Input : Mapped reads
# The  mapped  reads in bedfile format is summed  into  non-overlapping  10Kbp  bins  for  display  purposes.  
# Alternatively,  a  smoothing  filter is emploted,  comprised  of  a  moving  sum  with  window  size  of  10Kbp and a window slide of 100bp, followed by a moving median with window size of 10K bins and a window slide of a 100 bins.
def generate_PTR_dataframe(bedfile, outfile, OriC_coordinates, ter_coordinates):
    out_file = bedfile.replace(".bed", "_bins.csv")
    read_counts = []
    with open(bedfile) as fp:
        for line in fp:
            line = line.strip()
            line_split = line.split('\t')
            read_counts.append(int(line_split[2]))
    
    OriC_mean_depth =  statistics.mean(read_counts[(int(OriC_coordinates.split(',')[0]) - 1):(int(OriC_coordinates.split(',')[1]) - 1)])
    ter_mean_depth = statistics.mean(read_counts[(int(ter_coordinates.split(',')[0]) - 1):(int(ter_coordinates.split(',')[1]) - 1)])

    PTR_loci =  OriC_mean_depth / ter_mean_depth

    with open(outfile, 'w+') as out:
        out.write(str(bedfile.replace('.bed', ''))+' OriC_mean_depth :\t' + str(OriC_mean_depth) + '\n')
        out.write(str(bedfile.replace('.bed', ''))+' ter_mean_depthn:\t' + str(ter_mean_depth) + '\n')
        out.write(str(bedfile.replace('.bed', '')) + ' PTR_loci:\t' + str(PTR_loci) + '\n')
    out.close()
    

""" Start of Main Method/Pipeline """
if __name__ == '__main__':
    args = parser().parse_args()
    generate_PTR_dataframe(args.bedfile, args.outfile, args.OriC_coordinates, args.ter_coordinates)