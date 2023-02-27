#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:04:16 2023

@author: jy15
"""

import argparse
import csv
import os.path

import pandas as pd
from intervaltree import IntervalTree

# parse the input
parser = argparse.ArgumentParser(
    description="Alignment annotation and count genes, transcripts and biotypes"
)
parser.add_argument("--sam_file", type=str, required=True)
parser.add_argument("--outdir", type=str, required=True)
parser.add_argument(
    "--annotation_file",
    type=str,
    required=False,
    default="ensembl_transcripts_processed_primary_assembly.csv",
)
args = parser.parse_args()
if not (
    os.path.isfile(args.sam_file)
    and (args.sam_file[-4:] == ".sam" or args.sam_file[-4:] == ".SAM")
):
    print("Invalid SAM file input")
    exit()

if not os.path.isdir(args.outdir):
    print("Invalid output directory")
    exit()

# create a pandas dataframe with the processed gencode hg38 gff3
with open(args.annotation_file, newline="") as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
ensembl_transcripts = pd.read_csv(
    args.annotation_file,
    sep=dialect.delimiter,
    names=["Seqid", "Start", "End", "gene_name", "transcript_name", "transcript_type"],
)


# correct the gff3 seqid to match the RNAME in the SAM files with ensembl hg38
# reference release 107
def correct_gff3_Seqid_to_match_alignment(gff3_Seqid):
    if (gff3_Seqid[0:3]) == "chr":
        if gff3_Seqid[3] != "M":
            return gff3_Seqid[3:]
        else:
            return "MT"
    else:
        return gff3_Seqid


# function to construct an interval tree with pandas df containing features of one chr
# as the input
def construct_intervaltree(annotation_df_by_chr):
    temp_tree = IntervalTree()
    for i in range(annotation_df_by_chr.shape[0]):
        temp_line = annotation_df_by_chr.iloc[i]
        left = temp_line[1]
        right = temp_line[2]
        if left < right:
            temp_tree[left:right] = temp_line.values.tolist()
    return temp_tree


# construct a dictionary of interval trees with Seqid/chr as the key
transcript_IntervalTree_dict_by_chr = {}
all_Seqids = ensembl_transcripts.Seqid.unique().tolist()

for each_chr in all_Seqids:
    transcript_IntervalTree_dict_by_chr[
        correct_gff3_Seqid_to_match_alignment(each_chr)
    ] = construct_intervaltree(
        ensembl_transcripts.loc[ensembl_transcripts["Seqid"] == each_chr]
    )


# find the closest feature if there are overlapping features, set implementation
# return sequentially gene_name, transcript_name and transcript type, but not
# absolute_dist
def pickby_start_end_abs_distance_sum(query_start, query_end, feature_set):
    smallest_abs_dist = float("inf")
    for each in feature_set:
        temp_dist = abs(query_start - each[0]) + abs(query_end - each[1])
        if temp_dist < smallest_abs_dist:
            smallest_abs_dist = temp_dist
            smallest_abs_dist_match = each

    final_feature = smallest_abs_dist_match[2]
    return final_feature[-3], final_feature[-2], final_feature[-1]


# count with dictionary
def count_dict(count_key, count_dictionary):
    if count_key in count_dictionary:
        count_dictionary[count_key] += 1
    else:
        count_dictionary[count_key] = 1


# get all keys in the IntervalTree dictionary
all_Seqid_keys = set(transcript_IntervalTree_dict_by_chr.keys())

# count feature

gene_count: dict[str, int] = dict()  # count gene
transcript_count: dict[str, int] = dict()  # count transcript
biotype_count: dict[str, int] = dict()  # count biotype
unannotated_count = 0  # count unannotated aligned reads
unannotated_chr_read_count = (
    0  # count reads that are aligned to unannotated chr/reference
)
unannotated_chr: dict[
    str, int
] = dict()  # count alignment to a unannotated chr/reference

alignment_count = 0
with open(args.sam_file, "r") as sam_file:
    for alignment in sam_file:
        alignment_count += 1
        temp = alignment.strip().split("\t")
        legit_RNAME = temp[2]  # RNAME=temp[2]
        left = min(int(temp[3]), int(temp[7]))  # POS=int(temp[3]), PNEXT=int(temp[7])
        right = (
            left + abs(int(temp[8])) - 1
        )  # TLEN=int(temp[8]), abs(TLEN)=right-left+1
        if legit_RNAME in all_Seqid_keys:
            temp_match = transcript_IntervalTree_dict_by_chr[legit_RNAME][
                left:right
            ]  # tree search return a set of interval objects
            if len(temp_match) == 0:
                # add 1 to the unannotated count if no match found from tree search
                unannotated_count += 1
            else:
                (
                    match_gene,
                    match_transcript,
                    match_type,
                ) = pickby_start_end_abs_distance_sum(left, right, temp_match)
                count_dict(match_gene, gene_count)
                count_dict(match_transcript, transcript_count)
                count_dict(match_type, biotype_count)
        else:
            unannotated_count += 1
            unannotated_chr_read_count += 1
            count_dict(legit_RNAME, unannotated_chr)

biotype_count["unknow/unannotated"] = unannotated_count

# write the log file and the annotation and counting results
sample_name = args.sam_file.split(os.sep)[-1][:-4]

log_file = open(os.path.join(args.outdir, sample_name + "_log.txt"), "w")
log_file.writelines("Sample_name" + "\t" + sample_name + "\n")
log_file.writelines("Total_proper_read_pair_count" + "\t" + str(alignment_count) + "\n")
log_file.writelines(
    "Annotated_proper_read_pair_count"
    + "\t"
    + str(alignment_count - unannotated_count)
    + "\n"
)
log_file.writelines(
    "Unannotated_proper_read_pair_count" + "\t" + str(unannotated_count) + "\n"
)
log_file.writelines(
    "reads_mapped_to_ensembl_unannotated_reference"
    + "\t"
    + str(unannotated_chr_read_count)
    + "\n"
)
log_file.writelines("Total_gene_count" + "\t" + str(len(gene_count)) + "\n")
log_file.writelines("Total_transcript_count" + "\t" + str(len(transcript_count)) + "\n")
log_file.writelines(
    "Total_transcript_biotype_count" + "\t" + str(len(biotype_count)) + "\n"
)
log_file.close()


def write_files(file, input_list):
    f = open(file, "w")
    for each in input_list:
        temp = ""
        for i in range(len(each) - 1):
            temp += str(each[i])
            temp += "\t"
        temp += str(each[-1])
        f.writelines(temp + "\n")
    f.close()


sorted_gene_count_list = list(
    sorted(gene_count.items(), key=lambda item: item[1], reverse=True)
)
sorted_transcript_count_list = list(
    sorted(transcript_count.items(), key=lambda item: item[1], reverse=True)
)
sorted_biotype_count_list = list(
    sorted(biotype_count.items(), key=lambda item: item[1], reverse=True)
)
sorted_unannotated_chr_count_list = list(
    sorted(unannotated_chr.items(), key=lambda item: item[1], reverse=True)
)


write_files(
    os.path.join(args.outdir, sample_name + "_gene_count.txt"), sorted_gene_count_list
)
write_files(
    os.path.join(args.outdir, sample_name + "_transcript_count.txt"),
    sorted_transcript_count_list,
)
write_files(
    os.path.join(args.outdir, sample_name + "_biotype_count.txt"),
    sorted_biotype_count_list,
)
write_files(
    os.path.join(args.outdir, sample_name + "_unannotated_chr_count.txt"),
    sorted_unannotated_chr_count_list,
)
