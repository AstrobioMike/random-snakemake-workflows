#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding MAG-mapping workflow. The purpose is to generate MAG-level
coverage and detection information for each sample. It takes the contig-level coverage
and detection outputs output from bbmap's pileup.sh, and parses them to average all the values for 
all contigs in each MAG. 

E.g., if MAG_A has 2 contigs, and in Sample_X they look like this:

            coverage    detection
contig_1          20          0.8
contig_2          10          0.6

And MAG_B has 2 contigs, and in Sample_X they look like this:

            coverage    detection
contig_3           0          0.0
contig_4          10          0.2


The coverage output table would be like this:

        Sample_X
MAG_A         15
MAG_B          5

And the detection output table would be like this:

        Sample_X
MAG_A        0.7
MAG_B        0.1

## NOTE ##
It's written to work with a specific workflow right now, so there are some expectations about input filenames hard-coded in here.
##########

"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from Bio import SeqIO
from math import isnan
from numpy import NaN

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. The purpose is \
                                              generate MAG-level coverage and detection information for each sample. \
                                              See top of script for more details.")

required = parser.add_argument_group('required arguments')

required.add_argument("-s", "--sample-IDs-file", help="Single-column file holding unique sample IDs", action="store", required=True)
required.add_argument("-M", "--MAG-IDs-file", help="Single-column file holding unique MAG IDs", action="store", required=True)
required.add_argument("--MAG-directory", help="Directory holding MAG fastas", action="store", required=True)
required.add_argument("--sample-coverage-info-directory", help="Directory holding pileup.sh output tables", action="store", required=True)

parser.add_argument("--MAG-fasta-suffix", help='Suffix of MAG fasta files following unique MAG IDs (default: ".fa")', action="store", default=".fa")
parser.add_argument("-o", "--output-prefix", help='Desired output file prefix (default: "Combined-MAG-level")', action="store", default="Combined-MAG-level")
parser.add_argument("--output-directory", help='Directory for output files (default: "./")', action="store", default="./")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################

def main():

    MAG_IDs_list = file_to_list(args.MAG_IDs_file)
    sample_IDs_list = file_to_list(args.sample_IDs_file)

    header_dict = get_header_dict(MAG_IDs_list, args.MAG_fasta_suffix, args.MAG_directory)

    coverage_df, detection_df = build_cov_and_det_tabs(MAG_IDs_list, sample_IDs_list, header_dict, args.sample_coverage_info_directory)

    # writing out
    coverage_df_out_path = str(args.output_directory) + str(args.output_prefix) + "-coverages.tsv"
    coverage_df.to_csv(coverage_df_out_path, sep = "\t", index = False)

    detection_df_out_path = str(args.output_directory) + str(args.output_prefix) + "-detections.tsv"
    detection_df.to_csv(detection_df_out_path, sep = "\t", index = False)

################################################################################

### functions ###

def file_to_list(file):

    curr_list = [line.strip() for line in open(file)]

    return(curr_list)


def get_header_dict(MAG_IDs, fasta_suffix, MAG_directory):
    """ creates a dictionary of MAGs (keys) to their contig headers (values) """

    header_dict = {key: [] for key in MAG_IDs}

    for MAG in MAG_IDs:


        current_file = str(MAG_directory) + str(MAG) + str(fasta_suffix)

        with open(current_file) as in_fasta:

            for seq_record in SeqIO.parse(in_fasta, "fasta"):

                header_dict[MAG].append(seq_record.id)

    return(header_dict)


def read_in_sample_table(sample_ID, sample_coverage_info_directory):

    curr_file = str(sample_coverage_info_directory) + str(sample_ID) + "-pileup-contig-cov-and-det.tsv"

    curr_df = pd.read_csv(curr_file, sep = "\t", low_memory = False, header = 0, usecols = ["#ID", "Avg_fold", "Covered_percent"])

    # changing covered percent from percentages to proportions
    curr_df['Covered_percent'] = curr_df['Covered_percent'].apply(lambda x: x / 100)

    # renaming columns
    curr_df.rename(columns = {"#ID": "contig_ID", "Avg_fold": "cov", "Covered_percent": "det"}, inplace = True)

    return(curr_df)


def get_mean_cov_and_det_for_contigs(MAG_ID, sample_df, header_dict):
    """ gets the average coverage and detection for all contigs of a given MAG in a given sample """

    # subset table to just the contigs we want
    sub_df = sample_df[sample_df["contig_ID"].isin(header_dict[MAG_ID])]

    curr_cov_mean = sub_df["cov"].mean()
    curr_det_mean = sub_df["det"].mean()

    return(curr_cov_mean, curr_det_mean)


def build_cov_and_det_tabs(MAG_IDs, sample_IDs, header_dict, sample_coverage_info_directory):

    # initializing tables
    cov_df = pd.DataFrame(index = MAG_IDs)
    det_df = pd.DataFrame(index = MAG_IDs)

    # iterating through samples
    for sample in sample_IDs:

        curr_sample_df = read_in_sample_table(sample, sample_coverage_info_directory)

        # making dicts of coverages and detections of all MAGs in current sample
        curr_cov_dict = {}
        curr_det_dict = {}

        # iterating through MAGs in current sample
        for MAG in MAG_IDs:

            curr_cov_mean, curr_det_mean = get_mean_cov_and_det_for_contigs(MAG, curr_sample_df, header_dict)

            # adding values to current sample's dictionaries
            curr_cov_dict[MAG] = curr_cov_mean
            curr_det_dict[MAG] = curr_det_mean

        ## adding that sample's dictionaries to the building table
        cov_df[sample] = pd.Series(curr_cov_dict)
        det_df[sample] = pd.Series(curr_det_dict)

    ## moving index to first column and naming
    cov_df.reset_index(inplace = True)
    cov_df.rename(columns = {'index': 'MAG_ID'}, inplace = True)
    det_df.reset_index(inplace = True)
    det_df.rename(columns = {'index': 'MAG_ID'}, inplace = True)

    return(cov_df, det_df)


if __name__ == "__main__":
    main()
