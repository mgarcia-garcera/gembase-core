import argparse
from species_identifier import GembaseFunctions

def process_fasta(file_path, si: SpeciesIdentifier):
    output_path = file_path.rsplit('.', 1)[0] + ".renamed.fasta"

    with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
        current_header = ""
        current_sequence = []

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    # Write previous record
                    outfile.write(current_header + '\n')
                    for seq_line in current_sequence:
                        outfile.write(seq_line + '\n')

                current_sequence = []

                try:
                    parts = line[1:].split(maxsplit=1)
                    original_seqid = parts[0]
                    species_strain = parts[1] if len(parts) > 1 else ""

                    identifier = si.get_identifier_from_fasta_header(line)
                    current_header = f">{identifier} {species_strain} |OLD:{original_seqid}|"
                    print(f"{file_path}: {identifier}")
                except ValueError as e:
                    print(f"Warning in {file_path}: {e}")
                    current_header = line  # Use original if invalid

            else:
                current_sequence.append(line)

        # Write last entry
        if current_header:
            outfile.write(current_header + '\n')
            for seq_line in current_sequence:
                outfile.write(seq_line + '\n')

def main():
    parser = argparse.ArgumentParser(description="Generate species-strain identifiers.")
    parser.add_argument("fasta_files", nargs="+", help="FASTA file(s) to process")
    parser.add_argument("--history", required=True, help="Path to history TSV file")
    parser.add_argument("--reference", required=True, help="Path to species reference TSV file")

    args = parser.parse_args()

    si = SpeciesIdentifier(args.history, args.reference)

    for fasta_file in args.fasta_files:
        process_fasta(fasta_file, si)

if __name__ == "__main__":
    main()