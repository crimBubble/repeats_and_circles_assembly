import sys
from Bio import SeqIO

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py fasta_file output_file")
    sys.exit(1)

# Get the input fasta file and output file names from command-line arguments
fasta_file = sys.argv[1]
output_file = sys.argv[2]


# Function to extract the score from the read name
def extract_score(read_desc):
    info = read_desc.replace("(", "").replace(")", "")
    info = info.split("-")

    try:
        score = int(info[-1])
    except ValueError:
        score = float(read_desc.split("_")[-3])

    return score


# Read the fasta sequences and extract the scores
sequences = []
for record in SeqIO.parse(fasta_file, "fasta"):
    if len(record.seq) < 200:
        continue
    read_name = record.description
    score = extract_score(read_name)
    sequences.append((record, score))

# Sort the sequences based on the score in descending order
sorted_sequences = sorted(sequences, key=lambda x: x[1], reverse=True)

# Write the top 10 sequences with the highest scores to the output file
with open(output_file, "w") as output:
    for record, score in sorted_sequences[:10]:
        SeqIO.write(record, output, "fasta")
