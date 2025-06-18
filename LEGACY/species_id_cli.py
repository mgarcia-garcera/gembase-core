import argparse
from GembaseFunctions import SpeciesIdentifier

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
    parser = argparse.ArgumentParser(
        description="Assigns an identifier to a species assembly in a FASTA file."
    )
    parser.add_argument("fasta", help="Path to the input FASTA file (multi-FASTA allowed, same assembly/species expected)")
    parser.add_argument("-r", "--reference", default="species_reference.tsv", help="Path to reference file (species code)")
    parser.add_argument("-hi", "--history", default="species_history.tsv", help="Path to history file (species+strain identifiers)")
    parser.add_argument("-o", "--output", default="renamed_fastas", help="Output directory for renamed FASTA files")

    args = parser.parse_args()

    identifier = SpeciesIdentifier(
        reference_file=args.reference,
        history_file=args.history,
        output_dir=args.output
    )

    identifier.assign_identifier_to_fasta(args.fasta)

if __name__ == "__main__":
    main()