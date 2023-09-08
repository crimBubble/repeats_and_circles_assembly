#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO

# Note: Script has to be placed in the RepeatExplorer output directory!
in_file = "COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES_filter.csv"

# read COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES_filter.csv
# skip read number rows
# use SCL and CL information to create dict
dtype_dic = {'Cluster': str, 'Supercluster': str}
cluster_dict = pd.read_csv(in_file,
                           skiprows=0, index_col=0,
                           usecols=[0, 1], header=1,
                           dtype=dtype_dic,
                           delimiter='\t', ).squeeze().to_dict()

# get current working dir
current_directory = os.getcwd()
print(f"Current working dir: {current_directory}\n")

# step 1: copy SCL reads into new directory
print("Collect super cluster data...")

for cl_number in cluster_dict:

    # cl_number to string, add leading zeros
    dir_name = f"dir_CL{str(cl_number).zfill(4)}"

    print(f"Checking cluster: {dir_name}")

    # set location of clustering dir
    clustering_loc = os.path.join(current_directory, "seqclust", "clustering")

    # walk sub-folders in clustering and check for current CL
    # cur_dir = location, dirs = sub-folders, files = files
    for cur_dir, dirs, files in os.walk(clustering_loc):
        for directory in dirs:

            if directory == dir_name:
                # location of CL reads
                cl_reads_loc = os.path.join(cur_dir, directory, "reads.fasta")
                # target dir
                scl_dir_name = "SCL" + str(cluster_dict[cl_number]).zfill(3)
                # target location
                target_loc = os.path.join(current_directory, "ecc_assembly", scl_dir_name)
                # join target dir and location
                cl_reads_target = os.path.join(target_loc, f"{str(cl_number)}.reads.fasta")

                # create target dir if not present
                if not os.path.exists(target_loc):
                    print(f"Creating folder: {target_loc}")
                    os.makedirs(target_loc)

                # Copy reads
                # from "seqclust/clustering/clusters/dir_CLxxxx/reads.fasta"
                # to "/clust_assembly/SCLxxx/reads.fasta.x"
                print(f"Linking: {cl_reads_loc}")
                Path(cl_reads_target).symlink_to(Path(cl_reads_loc))

# step 2: merge reads from each SCL and create read list
print("Merge SCL data...")

# walk sub-folders in clust_assembly and merge reads
for cur_dir, dirs, files in os.walk(os.path.join(current_directory, "ecc_assembly")):
    if not len(files) > 0:
        # skip empty folders
        continue

    # get scl_dir_name from current dir
    scl_dir_name = os.path.basename(os.path.normpath(cur_dir))
    # strip SCL, keep number only
    scl_number = scl_dir_name[3:]

    # list with all reads
    all_records = []

    # read all fasta files in current dir
    for file in files:
        fasta_file = os.path.join(cur_dir, file)
        fasta_records = SeqIO.parse(fasta_file, "fasta")
        all_records.extend(fasta_records)
    # location and name of collected reads file
    collected_fasta = os.path.join(cur_dir, "reads.scl" + scl_number + ".fa")
    # write collected reads fasta file
    # print(f"Create collected reads file: {collected_fasta}")
    # SeqIO.write(all_records, collected_fasta, "fasta-2line")

    # collect all read names
    all_read_names = [record.id for record in all_records]

    # strip read orientation information (might be 1/2 or f/r)
    stripped_names = []
    last_read_name = []
    found_unknown_orientation = False
    for_ending = None
    rev_ending = None

    for name in all_read_names:
        if name.endswith("1") or name.endswith("2"):
            name = name[:-1]  # Remove the last character
            for_ending = 1
            rev_ending = 2
        elif name.endswith("for") or name.endswith("rev"):
            name = name[:-1]  # Remove the last character
            for_ending = "for"
            rev_ending = "rev"
        elif name.endswith("f") or name.endswith("r"):
            name = name[:-1]  # Remove the last character
            for_ending = "f"
            rev_ending = "r"
        else:
            last_read_name.append(name)
            found_unknown_orientation = True
            break

        stripped_names.append(name)

    if found_unknown_orientation:
        # found untypical orientation suffix, ask for user input
        print(f"Unknown read orientation suffix...")
        print(f"Currently processing read: {last_read_name}")
        for_ending_input = input(f"Enter suffix for forward orientation: ")
        rev_ending_input = input(f"Enter suffix for reverse orientation: ")
        print(f"Retry with new suffixes: for = {for_ending_input}; rev = {rev_ending_input}")

        for_ending = for_ending_input
        rev_ending = rev_ending_input

        stripped_names = []

        # re-run stripping of orientation information using user input
        for name in all_read_names:
            if name.endswith(for_ending_input) or name.endswith(rev_ending_input):
                suffix = len(for_ending_input)
                name = name[:-suffix]  # Remove the last characters

            stripped_names.append(name)

    # remove duplicates and keep only one occurrence of each read name
    unique_names = list(set(stripped_names))
    unique_names = sorted(unique_names)
    unique_names = [s[s.find('_')+1:s.rfind('_')] for s in unique_names]
    # for_ending = str(for_ending)
    # rev_ending = str(rev_ending)
    
    # remove prefix and suffix

    # file names for read lists
    collected_fasta_list = os.path.join(cur_dir, "reads.lst")
    collected_fasta_list_for = os.path.join(cur_dir, "reads.for.lst")
    collected_fasta_list_rev = os.path.join(cur_dir, "reads.rev.lst")

    # write read lists w/ and w/o read orientation info
    with open(collected_fasta_list, "w") as file:
        file.write("\n".join(unique_names))

    # with open(collected_fasta_list_for, "w") as file:
    #     file.write("\n".join([name + for_ending for name in unique_names]))

    # with open(collected_fasta_list_rev, "w") as file:
    #     file.write("\n".join([name + rev_ending for name in unique_names]))

    os.system(f"rm {cur_dir}/*.fasta")
