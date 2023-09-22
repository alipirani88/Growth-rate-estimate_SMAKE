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
    parser = argparse.ArgumentParser(description='Estimate ptr from coverage bed file')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-bedfile', action='store', dest="bedfile", help='Coverage Bed file', required=True)
    required.add_argument('-outfile', action='store', dest="outfile", help='Smoothed Coverage Depth', required=True)
    return parser

def generate_moving_sum_results(moving_sum_array, out_path):
    out_file = out_path + "_moving_sum_bins"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in moving_sum_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')
    peak = max(moving_sum_array)
    through = min([x for x in moving_sum_array if x !=0])
    PTR_moving = peak/through
    out_file_ptr = out_path.replace('.csv', '_PTR.txt')
    with open(out_file_ptr, 'w+') as out:
        out.write(str(out_path)+' moving_sum_array :\t'+str(PTR_moving)+'\n')

def calculate_per_of_reads(median_sliding_window_array, bedfile):
    stats_file = (bedfile.replace('.bed', '_alignment_stats')).replace("ptr", "picard")
    cmd = "grep \'mapped (\' " + stats_file + " | head -n1"
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    
    out = out.strip()
    out_split = str(out.decode('ascii')).split(' ')
    mapped_reads = out_split[0]
    out_file = os.path.dirname(bedfile) + "/" + os.path.basename(bedfile)[0:25] + "_perc_bins.csv"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                perc = (i*100)/int(mapped_reads)
                out.write(str(count)+','+str(perc)+'\n')
    return out_file

def generate_perc_coverage_graph(csv_matrix_file, ptr_value):
    out_path = csv_matrix_file.replace('.bed_perc_bins.csv', '_perc_coverage_graph.R')
    header = os.path.basename(csv_matrix_file).replace('.csv', '')
    with open(out_path, 'w') as out:
        print_string = "library(reshape)\nlibrary(ggplot2)\n"
        print_string = print_string + "dat=read.csv(\"%s\")\n" % (os.path.basename(csv_matrix_file))
        print_string = print_string + "mdf=melt(dat,id.vars=\"bin\")\n"
        print_string = print_string + "pdf(\"perc_coverage_graph.pdf\", width = 15, height = 10)\n"
        print_string = print_string + "p1 <- ggplot(mdf,aes(x=bin,y=value,colour=variable,group=variable)) + geom_line(size=1.5) + ggtitle(\"%s\")  + labs(x=\"Genomic Location\",y=\"Coverage in Percentage (perc reads/total)\") + theme(legend.title = element_blank(), legend.position=\"none\") + annotate(\"text\", x=max(mdf$bin)/2, y=max(mdf$value), label= \"PTR = %s\", size=5) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=10, face=\"bold\", margin = margin(10, 0, 10, 0)), axis.text.x = element_text(colour = \"black\", face= \"bold.italic\", vjust = 3)) + scale_color_manual(values=c(\"#3288bd\"))\n" % (header, ptr_value)
        print_string = print_string + "p1\ndev.off()"
        out.write(print_string)

def smoothing(read_counts, window, out_path, bedfile):
    print("Smoothing read depth data frame for PTR estimation")

    ## Step 1
    raw_count_sliding_window_array = []
    zero_bins = 0
    sixty_perc_bins = 0
    
    for i in range(0, len(read_counts), 100):
        start = i
        end = i + window
        if read_counts[start:end].count(0) == 10000:
            zero_bins += 1
        else:
            raw_count_sliding_window_array.append(read_counts[start:end])
        if read_counts[start:end].count(0) >= 6000:
                sixty_perc_bins += 1
   
    print("The number of bins with no mapped reads: %s" % str(zero_bins))
    
    ## Step 2
    moving_sum_array = []
    for i in raw_count_sliding_window_array:
        moving_sum_array.append(sum(i))
    generate_moving_sum_results(moving_sum_array, out_path)

    ## Step 3
    median_sliding_window_array = []
    for i in range(0, len(read_counts), 100):
        start = i
        end = i + window
        if len(moving_sum_array[start:end]) > 5000:
            median_sliding_window_array.append(statistics.median(moving_sum_array[start:end]))

    ## Step 4
    peak = max(median_sliding_window_array)
    through = min([x for x in median_sliding_window_array if x !=0])
    PTR_median = peak/through

    ## Step 5
    out_file_ptr = out_path.replace('.csv', '_PTR.txt')
    with open(out_file_ptr, 'a') as out:
        out.write(str(out_path)+' median_sliding_window_array :\t'+str(PTR_median)+'\n')
        out.write(str(out_path)+' Peak and trough location:\t' + str(median_sliding_window_array.index(peak)) + '\t' + str(median_sliding_window_array.index(through)) + '\n')
        out.write(str(out_path) + ' Peak and trough values:\t' + str(peak) + '\t' + str(through) + '\n')
    out.close()

    with open(out_path, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')

    # ## Step 6
    perc_bins_matrix = calculate_per_of_reads(median_sliding_window_array, bedfile)
    print(perc_bins_matrix)
    generate_perc_coverage_graph(perc_bins_matrix, PTR_median)

    return PTR_median



# Main Method
# The algorithm follows the procedure as described in this publication: Growth dynamics of gut microbiota in health and disease inferred from single metagenomic samples http://science.sciencemag.org/content/349/6252/1101.long with a few minor changes.
# Input : Mapped reads
# Steps: Input
# - The  mapped  reads in bedfile formatto  each  bacteria  was
# summed  into  non-
# overlap
# ping  10Kbp  bins  for  display  purposes.  Alternatively,  we
# employed  a  smoothing  filter,  comprised  of  a  moving  sum  with  window  size  of  10Kbp
# and a slide of 100bp, followed by a moving median with window size of 10K bins and a
# slide of a 100 bins.
def generate_PTR_dataframe(bedfile):
    out_file = bedfile.replace(".bed", "_bins.csv")
    window = 10000
    read_counts = []
    with open(bedfile) as fp:
        for line in fp:
            line = line.strip()
            line_split = line.split('\t')
            counts = int(line_split[2])
            read_counts.append(int(counts))
    
    PTR_median = smoothing(read_counts, window, out_file, bedfile)
    
    # perc_bins_matrix = os.path.dirname(bedfile) + "/" + os.path.basename(bedfile)[0:20] + "_perc_bins.csv"
    # generate_perc_coverage_graph(perc_bins_matrix, PTR_median, analysis)

""" Start of Main Method/Pipeline """
if __name__ == '__main__':
    args = parser().parse_args()
    generate_PTR_dataframe(args.bedfile)
