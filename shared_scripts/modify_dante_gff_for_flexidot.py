import sys

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py fasta_file output_file")
    sys.exit(1)

# Get the input fasta file and output file names from command-line arguments
gff_file = sys.argv[1]
output_file = sys.argv[2]


def process_gff(gff_file):
    with open(gff_file, 'r') as f:
        lines = f.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith("#"):
            modified_lines.append(line)  # Keep header lines intact
        else:
            fields = line.split("\t")
            attributes = fields[8].split(";")
            for i, attr in enumerate(attributes):
                if attr.startswith("Name="):
                    name_value = attr.split("=")[1]
                    # attributes[i] = name_value
                    fields[2] = name_value
                    break

            fields[8] = ";".join(attributes)
            modified_lines.append("\t".join(fields))

    modified_gff = "".join(modified_lines)
    return modified_gff


# Example usage
modified_gff = process_gff(gff_file)

with open(output_file, 'w') as f:
    f.write(modified_gff)

print("Modified GFF content has been written to domains.flexi.gff.")
