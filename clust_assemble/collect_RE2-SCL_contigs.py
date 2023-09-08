#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO

# Note: Script has to be placed in the RepeatExplorer output directory!
in_file = 'CLUSTER_TABLE_with_final_annotation.csv'

# read CLUSTER_TABLE_with_final_annotation.csv
# skip read number rows
# use SCL and CL information to create dict
dtype_dic = {'Cluster': str, 'Supercluster': str}
cluster_dict = pd.read_csv(in_file,
                           skiprows=6, index_col=0,
                           usecols=[0, 1], header=0,
                           dtype=dtype_dic,
                           delimiter='\t', ).squeeze().to_dict()

# get current working dir, should be inside RE2 default output file (same location as index.html)
current_directory = os.getcwd()
print("Current working directory: " + current_directory)

print("Create 'clust_assembly' folder and copying contig files...")

# loop over cluster directories in the cluster_dict and copy according to super clusters
for cluster_nr in cluster_dict:

    cl_dir_name = "dir_CL" + str(cluster_nr).zfill(4)
    clustering_dir = os.path.join(current_directory, "seqclust", "clustering")  # default RE2 naming scheme

    # loop over cluster directory and files
    for cur_dir, dirs, files in os.walk(clustering_dir):

        for directory in dirs:

            if directory == cl_dir_name:

                contigs_fa = os.path.join(cur_dir, directory, "contigs.info.fasta")

                if os.path.isfile(contigs_fa):

                    print(f"Found: {contigs_fa}")
                    # not all clusters contain contig files, copy tarean_contig instead

                else:

                    contigs_fa = os.path.join(cur_dir, directory, "tarean_contigs.fasta")

                    print(f"Found: {contigs_fa}")

                scl_dir_name = "SCL" + str(cluster_dict[cluster_nr]).zfill(3)
                target_dir = os.path.join(current_directory, "clust_assembly", scl_dir_name)
                target_location = os.path.join(target_dir, f"contigs.{str(cluster_nr)}.fasta")

                # create target directory
                if not os.path.exists(target_dir):
                    os.makedirs(target_dir)

                # linking contig files
                print(f"Linking: {contigs_fa}")
                Path(target_location).symlink_to(Path(contigs_fa))


print("Merging all contigs.X.fasta files into single fasta file...")
for cur_dir, dirs, files in os.walk(os.path.join(current_directory, "clust_assembly")):
    if not len(files) > 0:
        continue
    all_records = []

    for file in files:
        if file.startswith("contigs."):
            fasta_file = os.path.join(cur_dir, file)
            fasta_file_records = SeqIO.parse(fasta_file, "fasta")
            all_records.extend(fasta_file_records)

    fasta_location = os.path.join(cur_dir, "collected.contigs.fasta")

    print("Create file: " + fasta_location)
    SeqIO.write(all_records, fasta_location, "fasta-2line")
