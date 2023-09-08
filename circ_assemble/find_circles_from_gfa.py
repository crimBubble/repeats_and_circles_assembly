#!/usr/bin/env python

import argparse
import sys
import os

import pandas as pd
import networkx as nx
import re
import matplotlib.pyplot as plt

# Define global variable for input file
input_file_default = 'test_simple.gfa'


def parse_gfa_file(file):
    # read in input file and split by lines
    with open(file, 'r') as f:
        lines = f.read().splitlines()

    # create empty lists to store S and L rows
    s_rows = []
    l_rows = []

    # iterate over each line and check if it's an S or L row
    for line in lines:
        cols = line.split('\t')
        if cols[0] == 'S':
            # parse S row and append to s_rows list
            seq_length = int(re.sub('(?i)ln:i:', '', cols[3]))
            seq_coverage = float(re.sub('(?i)dp:f:', '', cols[4]))

            row = {'type': cols[0],
                   'seqid': (cols[1]),
                   'seq': cols[2],
                   'len': seq_length,
                   'cov': seq_coverage}

            s_rows.append(row)

        elif cols[0] == 'L':
            # parse L row and append to l_rows list
            row = {'type': cols[0],
                   'from_id': cols[1],
                   'from_strand': cols[2],
                   'to_id': cols[3],
                   'to_strand': cols[4],
                   'distance': str(cols[5])}
            l_rows.append(row)

    # create pandas dataframes for S and L rows
    s_df = pd.DataFrame(s_rows, columns=['type', 'seqid', 'seq', 'len', 'cov'])
    l_df = pd.DataFrame(l_rows, columns=['type', 'from_id', 'from_strand', 'to_id', 'to_strand', 'distance'])

    return s_df, l_df


def filter_gfa_cov(df_seq, df_link, min_cov: float):
    # filter Sequences df based on cov threshold
    filtered_seqs = df_seq[df_seq['cov'] >= min_cov]

    # filter Linkages df to only contain from_id and to_id values that are still present in the seqid column of the
    # Sequences df
    filtered_links = df_link[
        df_link['from_id'].isin(filtered_seqs['seqid']) & df_link['to_id'].isin(filtered_seqs['seqid'])]

    return filtered_seqs, filtered_links


def filter_gfa_deadends(df_seq_cov, df_link_cov):
    # filter dead ends (contigs with nothing linked to one contig-end)
    # detected small circles of just a single contig

    dead_end_ids = []
    more_changes = True

    for seqid in df_seq_cov['seqid']:

        # display output meaning "hi, I am still running"
        for i in range(4):
            print('\r' + ' ' * (3 - i) + '...' + ' ' * i + f' run filtering on {seqid}', end='')

        linkages = df_link_cov[(df_link_cov['from_id'] == seqid) | (df_link_cov['to_id'] == seqid)]

        not_deadend = True

        while not_deadend:

            if len(linkages) == 0:  # no linkages, dead ends
                dead_end_ids.append(seqid)
                not_deadend = False

            elif len(linkages) == 1:  # one linkage, potential small circle
                if linkages['from_id'].iloc[0] == linkages['to_id'].iloc[0]:  # small circle, detect again later
                    break

                else:
                    dead_end_ids.append(seqid)  # dead ends, linkage to only one other contig
                    not_deadend = False

            else:
                left_end = False
                right_end = False

                for i, linkage in linkages.iterrows():
                    if linkage['from_id'] == seqid and linkage['from_strand'] == '+':
                        left_end = True
                    elif linkage['to_id'] == seqid and linkage['to_strand'] == '-':
                        left_end = True
                    elif linkage['from_id'] == seqid and linkage['from_strand'] == '-':
                        right_end = True
                    elif linkage['to_id'] == seqid and linkage['to_strand'] == '+':
                        right_end = True

                if left_end and right_end:
                    break

                else:
                    dead_end_ids.append(seqid)
                    not_deadend = False

    filtered_seqs = df_seq_cov[~df_seq_cov['seqid'].isin(dead_end_ids)]

    filtered_links = df_link_cov[
        df_link_cov['from_id'].isin(filtered_seqs['seqid']) & df_link_cov['to_id'].isin(filtered_seqs['seqid'])]

    if filtered_seqs.equals(df_seq_cov):
        more_changes = False

    print('\n')

    return filtered_seqs, filtered_links, more_changes


def detect_small_circles(df_seq_filter, df_link_filter, assembly_graph):
    # Detect circles containing exactly two contigs

    small_circles_one = []
    small_circles_two = []

    circ_patterns = [['L', 'id1', '-', 'id2', '+', '0M'],
                     ['L', 'id1', '+', 'id2', '-', '0M'],
                     ['L', 'id1', '+', 'id2', '+', '0M'],
                     ['L', 'id1', '-', 'id2', '-', '0M'],
                     ['L', 'id1', '+', 'id2', '-', '0M'],
                     ['L', 'id1', '-', 'id2', '+', '0M'],
                     ['L', 'id1', '+', 'id2', '+', '0M'],
                     ['L', 'id2', '+', 'id1', '+', '0M'],
                     ['L', 'id1', '-', 'id2', '+', '0M'],
                     ['L', 'id2', '+', 'id1', '-', '0M'],
                     ['L', 'id1', '+', 'id2', '-', '0M'],
                     ['L', 'id2', '-', 'id1', '+', '0M']]

    L_header = ['type', 'from_id', 'from_strand', 'to_id', 'to_strand', 'distance']

    pd_small_circles_two_patterns = pd.DataFrame(circ_patterns, columns=L_header)

    for seqid in df_seq_filter['seqid']:

        # display output meaning "hi, I am still running"
        for i in range(4):
            print('\r' + ' ' * (3 - i) + '...' + ' ' * i + f' detect small circles including contig {seqid}', end='')

        linkages = df_link_filter[(df_link_filter['from_id'] == seqid) | (df_link_filter['to_id'] == seqid)]

        # detect single contig circles
        if len(linkages) == 1:

            if linkages['from_id'].iloc[0] == linkages['to_id'].iloc[0]:

                print(f'\nCircle from single contig found.')

                small_circles_one.append([int(seqid)])

        # detect two contig circles
        else:
            adj_seqids = list(assembly_graph.neighbors(seqid))

            for adj_seqid in adj_seqids:

                adj_links = linkages[(linkages['from_id'] == adj_seqid) | (linkages['to_id'] == adj_seqid)]

                expected_patterns = pd_small_circles_two_patterns.replace({'id1': seqid, 'id2': adj_seqid})

                if len(adj_links) >= 2:

                    overlaps = adj_links.merge(expected_patterns, how='left', indicator=True)
                    overlaps_found = (overlaps['_merge'] == 'both').all()

                    if overlaps_found:
                        print(f'\nSmall circle with two contigs found.')
                        small_circles_two.append([seqid, adj_seqid])

    print('\n')

    if small_circles_two:
        print(f'Removing duplicates from small circles.')
        # make it tuples and filter unique pairs (order does not matter)
        small_circles_two = set(tuple(sorted(set(pair))) for pair in small_circles_two)
        # transform back to list of lists
        small_circles_two = [list(t) for t in small_circles_two]

    return small_circles_one, small_circles_two


def create_assembly_graph(df_seq, df_link):
    # Create a new empty graph
    assembly_graph = nx.Graph()

    # Add nodes to the graph
    nodes = [(row['seqid']) for idx, row in df_seq.iterrows()]
    assembly_graph.add_nodes_from(nodes)

    # Add edges to the graph
    edges = [(row['from_id'], row['to_id']) for idx, row in df_link.iterrows()]
    assembly_graph.add_edges_from(edges)

    return assembly_graph


def filter_circles(circles, contig):
    return [filtered for filtered in circles if contig in filtered]


def export_circ_seq(circles, df_seq_filter, in_file):
    # collect sequences from circle lists and export as fasta

    output_file = in_file.replace('.gfa', '.circles.fasta')
    print(f'Output file: {output_file}')

    with open(output_file, 'w') as out_file:

        for circle in circles:
            # make a subset of df_seq to only contain seqids from circle
            # seqid in circle is int and str in df_seq, convert
            select_df_seq = df_seq_filter.loc[df_seq_filter['seqid'].isin([str(i) for i in circle])]
            seq_seq = ''.join(select_df_seq['seq'].tolist())

            # weighted coverage
            seq_len = select_df_seq['len'].sum()
            seq_cov = (select_df_seq['len'] * select_df_seq['cov']).sum() / seq_len

            # seq name
            seq_name = f'>circle_nodes_' + ':'.join(str(node) for node in circle) + \
                       f' length={seq_len}' + \
                       f' depth={seq_cov:.2f}x' + \
                       ' circular=true'

            out_file.write(seq_name + os.linesep)
            out_file.write(seq_seq + os.linesep)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Parse GFA file and receive circular patterns')
    parser.add_argument('--input', '-i', type=str, default=input_file_default,
                        help='Input GFA file (default: test_simple.gfa)')
    parser.add_argument('--coverage', '-c', type=float, default=0.01,
                        help='Minimal coverage depth of a contig. Contigs with lesser coverage will be filtered out ('
                             'default: 0.01)')
    parser.add_argument('--display', action='store_true', help='Display gfa network. WARNING: Circle detection will '
                                                               'not continue while the contig network is displayed.')
    args = parser.parse_args()

    print('Running with:', ' '.join(sys.argv[1:]))

    input_file = args.input

    try:
        input_file_size = os.path.getsize(input_file)
    except FileNotFoundError:
        print(f'Input file {input_file} not found.')
        input_file_size = 0

    if input_file_size == 0:
        input_file_new = os.path.join(os.path.dirname(input_file), '004_bridges_applied.gfa')
        try:
            input_file_size = os.path.getsize(input_file_new)
            input_file = input_file_new
        except FileNotFoundError:
            print(f'Input alternative file {input_file_new} not found.')
            input_file_size = 0

    if input_file_size == 0:
        input_file_new = os.path.join(os.path.dirname(input_file), '003_overlaps_removed.gfa')
        try:
            os.path.getsize(input_file_new)
            input_file = input_file_new
        except FileNotFoundError:
            print(f'Input alternative file {input_file_new} not found.')
            exit(1)

    print(f'Reading and filtering {input_file} ...')

    # Parse GFA file
    df_seq, df_link = parse_gfa_file(input_file)

    df_seq_filter, df_link_filter = filter_gfa_cov(df_seq, df_link, min_cov=args.coverage)

    more_changes = True

    while more_changes:
        df_seq_filter, df_link_filter, more_changes = filter_gfa_deadends(df_seq_filter, df_link_filter)

    # re-run parsing and filtering with bridges file if nothing is left after filtering
    if df_seq_filter.empty:
        input_file = os.path.join(os.path.dirname(input_file), '004_bridges_applied.gfa')

        print(f'Reading and filtering {input_file} ...')

        # Parse GFA file
        df_seq, df_link = parse_gfa_file(input_file)

        df_seq_filter, df_link_filter = filter_gfa_cov(df_seq, df_link, min_cov=args.coverage)

        more_changes = True

        while more_changes:
            df_seq_filter, df_link_filter, more_changes = filter_gfa_deadends(df_seq_filter, df_link_filter)

    # re-run parsing and filtering with overlaps file if nothing is left after filtering
    if df_seq_filter.empty:
        input_file = os.path.join(os.path.dirname(input_file), '003_overlaps_removed.gfa')

        print(f'Reading and filtering {input_file} ...')

        # Parse GFA file
        df_seq, df_link = parse_gfa_file(input_file)

        df_seq_filter, df_link_filter = filter_gfa_cov(df_seq, df_link, min_cov=args.coverage)

        more_changes = True

        while more_changes:
            df_seq_filter, df_link_filter, more_changes = filter_gfa_deadends(df_seq_filter, df_link_filter)

    print(f'Sequences (head): \n{df_seq_filter.head()}')
    print(f'Linkages (head): \n{df_link_filter.head()}')

    # Write filtered dataframes to output file
    output_file = input_file.replace('.gfa', '.filter.gfa')
    print(f'Saving filtered GFA file: {output_file}')
    with open(output_file, 'w') as f:
        for idx, row in df_seq_filter.iterrows():
            f.write(f"{row['type']}\t{row['seqid']}\t{row['seq']}\tln:i:{row['len']}\tdp:f:{row['cov']}\n")
        for idx, row in df_link_filter.iterrows():
            f.write(f"{row['type']}\t{row['from_id']}\t{row['from_strand']}\t{row['to_id']}\t{row['to_strand']}\t{row['distance']}\n")

    print(f'Creating contigs network ...')
    assembly_graph = create_assembly_graph(df_seq_filter, df_link_filter)

    print(f'Starting circle detection ...')

    print(f'Detecting small circles ...')
    small_circles_one, small_circles_two = detect_small_circles(df_seq_filter, df_link_filter, assembly_graph)
    print(f'Small circles (single contig): {small_circles_one}')
    print(f'Small circles (two contigs): {small_circles_two}')

    if args.display:
        nx.draw(assembly_graph, with_labels=True)
        plt.show()

    print(f'Detecting circles ...')
    best_contig = list(assembly_graph.nodes)[0]  # best assembled contig is usually contig 1
    long_circles = nx.cycle_basis(assembly_graph, root=best_contig)

    # filter for circles to contain contig '1'
    long_circles_filter = filter_circles(long_circles, '1')
    # transform lists from strings to integer to use in bandage
    long_circles_filter = [[int(x) for x in sublist] for sublist in long_circles_filter]

    print('Circular patterns:')

    if len(long_circles_filter) == 0:
        print(f'No large circles found.')

    else:
        for circ in long_circles_filter:
            print(circ)

    print('Collecting circles ...')

    circles = []

    if small_circles_one:
        circles += small_circles_one

    if small_circles_two:
        circles += small_circles_two

    if long_circles_filter:
        circles += long_circles_filter

    clean_circles = []
    for sub_circles in circles:
        if sub_circles not in clean_circles:
            clean_circles.append(sub_circles)

    print('Exporting circular sequences ...')

    export_circ_seq(clean_circles, df_seq_filter, input_file)


if __name__ == '__main__':
    main()
