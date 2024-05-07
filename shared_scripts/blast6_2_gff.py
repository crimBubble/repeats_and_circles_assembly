#!/usr/bin/env python
import pandas as pd
from typing import Union

default_m6_header = "6 delim=\\t qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"


def check_m6_header(infmt):
    """
    Check m6 header
    delimiter ~ separator
    :option 1: default header (-outfmt 6)
    :option 2: default header with different separator
    :option 3: custom header with default separator, minimum requirements for proper gff are id start end
    :option 4: custom header with custom separator
    """

    while True:

        # read header, first element should be always 6
        m6_header = infmt.split(" ")

        # default header with no user input or user input equal to default
        if m6_header == default_m6_header.split(" "):
            m6_sep = m6_header[1].strip("delim=")
            m6_header = m6_header[2:]

            print('Using default input format.')

            break

        # header does not start with 6, results in error
        if m6_header[0] != "6":
            print(f'Unexpected input format: {m6_header[0]}\n'
                  f'Expected "-i/-infmt" starts with "6 ..."')
            exit(1)

        # header does start with 6, but is only 6, ask for to continue with default settings
        if m6_header[0] == "6" and len(m6_header) == 1:

            user_answer = input(f'-infmt "6" does not need to be specified. '
                                f'Continue with default input format?')

            if user_answer in ["y", "yes", "", "Yes", "Y", "YES"]:

                m6_header = default_m6_header.split(" ")[2:]
                m6_sep = default_m6_header.split(" ")[1].strip("delim=")

                break

            else:
                exit(0)

        # header does start with 6, remove 6 and continue checks
        else:
            m6_header = m6_header[1:]

        # short-1 header only includes custom delimiter (or default delimiter)
        if len(m6_header) == 1:

            if m6_header[0].startswith("delim="):

                m6_sep = m6_header[0].strip("delim=")
                m6_header = default_m6_header.split(" ")[2:]

                print(f"Using default blastn header with delimiter {m6_sep}")

                break

            else:
                print("Improper values in input format given.")
                exit(1)

        # header does not contain enough values for a proper gff output
        if len(m6_header) < 3:
            print('Not enough values in input format given.\n'
                  'Provide at least "6 sseqid sstart ssend" or "6 qseqid qstart qend".')
            exit(1)

        # function to check if all 3 minimum requirements are in custom header
        def check_min_req_val(header):

            id_true, start_true, end_true = False, False, False

            if any("id" in s for s in header):
                id_true = True

            if any("start" in s for s in header):
                start_true = True

            if any("end" in s for s in header):
                end_true = True

            if id_true and start_true and end_true:
                return True

        # custom header with default delimiter
        if len(m6_header) >= 3 and not m6_header[0].startswith("delim="):

            if check_min_req_val(m6_header):

                m6_sep = default_m6_header.split(" ")[1].strip("delim=")

                break

            else:
                print('Improper values in input format given.\n'
                      'Provide at least "6 sseqid sstart ssend" or "6 qseqid qstart qend".')
                exit(1)

        # custom header with custom delimiter
        if len(m6_header) >= 3 and m6_header[0].startswith("delim="):

            if check_min_req_val(m6_header):

                m6_sep = m6_header[0].strip("delim=")
                m6_header = m6_header[1:]

                break

            else:
                print('Improper values in input format given.\n'
                      'Provide at least "6 delim= sseqid sstart ssend" or "6 delim= qseqid qstart qend".')
                exit(1)

    return m6_header, m6_sep


def check_mode_type(gff_mode, gff_type, m6_header):
    # set is_mode_q to true or false, and set default type if not set by user
    if gff_mode == "q":
        is_mode_q = True
        if gff_type is None:
            gff_type = "sseqid"

    elif gff_mode == "s":
        is_mode_q = False
        if gff_type is None:
            gff_type = "qseqid"

    else:
        print(f'Unknown mode. Please use --mode q or --mode s')
        exit(1)

    # depending on mode, check type and if correct start and end values are given
    if is_mode_q:

        if gff_type == "qseqid":
            print(f'--mode {gff_mode} does not work with --type qseqid')
            exit(1)

        if "sseqid" not in m6_header and gff_type == "sseqid":
            print(f'--mode {gff_mode} does need "sseqid" to work')
            exit(1)

        if "sstart" not in m6_header:
            print(f'--mode {gff_mode} does need "sstart" to work')
            exit(1)

        if "send" not in m6_header:
            print(f'--mode {gff_mode} does need "send" to work')
            exit(1)

        else:
            pass

    if not is_mode_q:

        if gff_type == "sseqid":
            print(f'--mode {gff_mode} does not work with --type sseqid')
            exit(1)

        if "qseqid" not in m6_header and gff_type == "qseqid":
            print(f'--mode {gff_mode} does need "qseqid" to work')
            exit(1)

        if "qstart" not in m6_header:
            print(f'--mode {gff_mode} does need "qstart" to work')
            exit(1)

        if "qend" not in m6_header:
            print(f'--mode {gff_mode} does need "qend" to work')
            exit(1)

        else:
            pass

    return gff_mode, gff_type


def read_m6(m6_file, m6_header, m6_sep):
    """
    read the m6 file
    :input 1: m6 file from blast(n) with -outfmt 6
    :input 2: minimum bit score defined by user for filtering
    :input 3: minimum e-value defined by user for filtering
    :return: blast hits in gff format
    """

    # read blast file
    print(f'Reading: {m6_file}')

    blast_in = pd.read_csv(m6_file, index_col=None, header=None, delimiter=m6_sep, engine='python')

    # check if header length is equal to number of columns
    col_number = len(blast_in.iloc[1])
    header_number = len(m6_header)

    if col_number != header_number:
        print(f'Number of elements in blast header: {header_number} '
              f'does not match number of columns in file: {col_number}'
              f'\nCheck for double spaces.')
        exit(1)

    else:
        blast_in.columns = m6_header

    print(f'Imported...\n'
          f'{blast_in.head(3)}\n'
          f'...\n'
          f'{blast_in.tail(3).to_string(header=False)}')

    return blast_in


def filter_m6(blast_input: pd.DataFrame,
              maximum_e_value: Union[float, None],
              minimum_bit_score:  Union[float, None]) -> pd.DataFrame:

    # Function to filter for minimum e-value or bitscore if they are not None

    if maximum_e_value is not None:
        blast_input_filtered = blast_input.loc[blast_input['evalue'] <= maximum_e_value]
    else:
        blast_input_filtered = blast_input

    if minimum_bit_score is not None:
        blast_input_filtered = blast_input_filtered.loc[blast_input_filtered['bitscore'] >= minimum_bit_score]
    else:
        blast_input_filtered = blast_input_filtered

    return blast_input_filtered


def convert_m6(blast_input, gff_mode, gff_type, m6_header):
    if gff_mode == "q":
        is_mode_q = True
    else:
        is_mode_q = False

    # build raw gff column by column

    # seqid, initialize new data frame
    if is_mode_q:
        gff_raw = blast_input[["qseqid"]].copy()
    else:
        gff_raw = blast_input[["sseqid"]].copy()

    gff_raw.rename(columns={gff_raw.columns[0]: 'seqid'}, inplace=True)

    # source
    gff_raw["source"] = "blast"

    # type
    if is_mode_q:
        if gff_type != "sseqid":
            blast_input[gff_type] = gff_type

        gff_raw["type"] = blast_input[gff_type]

    else:
        if gff_type != "qseqid":
            blast_input[gff_type] = gff_type

        gff_raw["type"] = blast_input[gff_type]

    # start
    if is_mode_q:
        gff_raw["start"] = blast_input["qstart"]
    else:
        gff_raw["start"] = blast_input["sstart"]

    # end
    if is_mode_q:
        gff_raw["end"] = blast_input["qend"]
    else:
        gff_raw["end"] = blast_input["send"]

    # score
    gff_raw["score"] = "."
    # gff_raw["score"] = blast_input["bitscore"]

    # strand
    if gff_raw["start"] <= gff_raw["end"]:
        gff_raw["strand"] = "+"
    else:
        # If 'start' is greater than 'end', swap them and set the strand to '-'
        temp_start = gff_raw["start"]
        gff_raw["start"] = gff_raw["end"]
        gff_raw["end"] = temp_start
        gff_raw["strand"] = "-"

    # phase
    gff_raw["phase"] = "."

    # attribute
    def combine_columns(row, col):
        new_line = str()
        for i in range(0, len(col) - 1):
            new_line += f'{col[i]}={str(row[col[i]])};'
        return new_line

    gff_raw["attributes"] = blast_input.apply(combine_columns, axis=1, col=m6_header)

    return gff_raw


def write_gff(m6_file, gff_all_hits_list_sorted):
    """
      """

    if m6_file.endswith(".m6"):
        gff_file = m6_file.rstrip(".m6") + ".gff"
    else:
        gff_file = m6_file + ".gff"

    print(f'Gff file is: {gff_file}')

    handle_gff = open(gff_file, "w")
    handle_gff.write("##gff-version 3\n")
    handle_gff.write("##this gff works with the original genome fasta\n\n")

    gff_all_hits_list_sorted.to_csv(path_or_buf=handle_gff, sep='\t', na_rep='', float_format=None, columns=None,
                                    header=False, index=False, index_label=None, mode='a', encoding=None, decimal='.')

    handle_gff.close()


def main(args):
    # get user input
    m6_file = args.file
    minimum_bit_score = args.min_bitscore
    maximum_e_value = args.max_evalue

    # check --infmt
    print('Checking user input...')

    m6_header, m6_sep = check_m6_header(args.infmt)

    print(f'Blast header is: {m6_header}\n'
          f'Column delimiter is: {m6_sep}')

    # check --mode and --type
    gff_mode, gff_type = check_mode_type(args.mode, args.type, m6_header)
    print(f'Mode is: {gff_mode}\n'
          f'Feature type is: {gff_type}\n')

    blast_input = read_m6(m6_file, m6_header, m6_sep)

    # filter for minimum scores
    print(f'Filtering blast table...\n'
          f'Minimum e-values is: {maximum_e_value}\n'
          f'Minimum bitscore is: {minimum_bit_score}\n')
    blast_input_fltr = filter_m6(blast_input, maximum_e_value, minimum_bit_score)

    # Converting blast table to GFF
    print(f'Converting blast table to gff3 format')
    gff_df_fltr_srt = convert_m6(blast_input_fltr, gff_mode, gff_type, m6_header)

    # Writing new gff file
    print(f'Writing new gff file...')
    write_gff(m6_file, gff_df_fltr_srt)

    print(f'Done.')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert a tabular (6) blast file to gff3 compatible format."
    )

    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument("-f",
                               "--file",
                               type=str,
                               required=True,
                               help="Blast file to convert.")

    parser.add_argument("-i",
                        "--infmt",
                        type=str,
                        default=default_m6_header,
                        help=f"Input format equivalent to -outfmt 6 of blastn.")

    parser.add_argument("-m",
                        "--mode",
                        type=str,
                        default="q",
                        help=f"Create GFF file to annotate query (q) or subject/db (s).")

    parser.add_argument("-t",
                        "--type",
                        type=str,
                        default=None,
                        help=f'Type of feature in GFF output file (see GFF3 file format). '
                             f'Dependent on --mode (sseqid for q, qseqid for s). '
                             f'Or use custom e.g., "feature" or "repeat"')

    parser.add_argument("-s",
                        "--min_bitscore",
                        type=float,
                        default=None,
                        help="Minimum bitscore for filtering blast hits. "
                             "Bitscore needs to be existent in blast file.")

    parser.add_argument("-e",
                        "--max_evalue",
                        type=float,
                        default=None,
                        help="Maximum e-value for filtering blast hits. "
                             "E-value needs to be existent in blast file.")

    args = parser.parse_args()
    main(args)
