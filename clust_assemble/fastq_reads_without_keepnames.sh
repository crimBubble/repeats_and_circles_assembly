#!/bin/bash

# Store the loop variable from the command-line argument, and original read files
dir="$1"
READ_1="$2"
READ_2="$3"

# Perform operations using the loop variable
echo "Processing directory: ${dir}"

# Extract sequence information from the fasta file
grep -v '^>' "${dir}reads.for.fasta" > "${dir}reads.for.seq"
grep -v '^>' "${dir}reads.rev.fasta" > "${dir}reads.rev.seq"

# Search for sequences in original fastq files using the sequences instead of the read names
# -m 0 to allow for no mismatches, still this will harvest duplicates (removed later)
seqkit grep -s -f "${dir}reads.for.seq" \
 -o "${dir}init.reads.for.fastq" -w 0 -j 1 -m 0 -P "${READ_1}"
seqkit grep -s -f "${dir}reads.rev.seq" \
 -o "${dir}init.reads.rev.fastq" -w 0 -j 1 -m 0 -P "${READ_2}"

# Make a read name list to check for pairs
# grep "^@" will also take quality information, should be removed in next step
# remove starting @ to make real list, remove description (might contain orientation which leads to no output)
grep "^@" "${dir}init.reads.for.fastq" | sed 's/^@//' | sed 's/ .*//' > "${dir}tmp1.txt"
grep "^@" "${dir}init.reads.rev.fastq" | sed 's/^@//' | sed 's/ .*//' > "${dir}tmp2.txt"

# Compare read name lists
grep -xF -f "${dir}tmp1.txt" "${dir}tmp2.txt" > "${dir}pe1.reads.fastq.lst"

# Collect paired-end reads
seqtk subseq "${dir}init.reads.for.fastq" "${dir}pe1.reads.fastq.lst" > \
 "${dir}pe1.reads.for.fastq"
seqtk subseq "${dir}init.reads.rev.fastq" "${dir}pe1.reads.fastq.lst" > \
 "${dir}pe1.reads.rev.fastq"

# Remove duplicates
seqkit rmdup "${dir}pe1.reads.for.fastq" -s -P -o "${dir}rmdup.reads.for.fastq"
seqkit rmdup "${dir}pe1.reads.rev.fastq" -s -P -o "${dir}rmdup.reads.rev.fastq"

# Check again for paired-end reads
grep "^@" "${dir}rmdup.reads.for.fastq" | sed 's/^@//' | sed 's/ .*//' > "${dir}tmp1.txt"
grep "^@" "${dir}rmdup.reads.rev.fastq" | sed 's/^@//' | sed 's/ .*//' > "${dir}tmp2.txt"
grep -xF -f "${dir}tmp1.txt" "${dir}tmp2.txt" > "${dir}pe2.reads.fastq.lst"

# Collect final PE reads
seqtk subseq "${dir}init.reads.for.fastq" "${dir}pe2.reads.fastq.lst" > \
"${dir}reads.for.fastq"
seqtk subseq "${dir}init.reads.rev.fastq" "${dir}pe2.reads.fastq.lst" > \
"${dir}reads.rev.fastq"

# Remove temporary files
rm "${dir}tmp1.txt" "${dir}tmp2.txt"
rm "${dir}pe1.reads.fastq.lst" "${dir}pe2.reads.fastq.lst"
rm "${dir}pe1.reads.for.fastq" "${dir}pe1.reads.rev.fastq"
rm "${dir}init.reads.for.fastq" "${dir}init.reads.rev.fastq"
rm "${dir}rmdup.reads.for.fastq" "${dir}rmdup.reads.rev.fastq"

echo "Done processing directory: ${dir}"
