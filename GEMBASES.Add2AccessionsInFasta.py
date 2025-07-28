import re
import argparse

parser = argparse.ArgumentParser(description="Cluster species FASTAs with MMseqs2")
parser.add_argument("-i", "--input", required=True, help="Species identifier (e.g., AABB.XXX)")
parser.add_argument("-o", "--output", required=True, help="Path to folder with FASTA files")
parser.add_argument("-a", "--add", required=False, default="mmseqs", help="Clustering method")

args = parser.parse_args()

accession_pattern = re.compile(r'⁽>([\w.]+)')

input_file = args.input
output_file = args.output

what2add = args.add

if input_file == output_file:
    print(f"input file should not be the same as the output")
else:
    with open(input_file, 'r') as infile, open(output_file, 'a') as outfile:
        for line in infile:
            if line.startswith('>'):
                match = accession_pattern.match(line)
                if match:
                    accession = match.group(1)
                    modified = accession + what2add + line[len(accession):]
                    outfile.write(modified)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

