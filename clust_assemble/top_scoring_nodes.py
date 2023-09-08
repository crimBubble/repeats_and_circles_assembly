import sys
from Bio import SeqIO

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py fasta_file output_file")
    sys.exit(1)

# Get the input fasta file and output file names from command-line arguments
fasta_file = sys.argv[1]
output_file = sys.argv[2]


# Function to calculate the score for a sequence
def calculate_score(seq):
    length = len(seq)
    coverage = float(seq.description.split("_")[-1])
    return length * coverage


# Read the fasta sequences and calculate the scores
sequences = []
for record in SeqIO.parse(fasta_file, "fasta"):
    if len(record.seq) < 200:
        continue
    score = calculate_score(record)
    sequences.append((record, score))

# Sort the sequences based on the score in descending order
sorted_sequences = sorted(sequences, key=lambda x: x[1], reverse=True)

# Write the top 10 sequences with the highest scores to the output file
with open(output_file, "w") as output:
    for record, score in sorted_sequences[:10]:
        SeqIO.write(record, output, "fasta")
