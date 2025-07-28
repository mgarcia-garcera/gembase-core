import sys
import os
import tempfile

def parse_taxid_mapping(mapping_file):
    taxid_map = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                species_id = parts[0]
                taxids = parts[1].split()
                taxid_map[species_id] = taxids[0]  # Use only the first taxid
    return taxid_map

def update_fasta_headers(fasta_file, tmp_file, taxid_map):
    # Write to a temporary file first
    with open(tmp_file, 'w') as tmp_out:
        with open(fasta_file, 'r') as fin:
            for line in fin:
                if line.startswith('>'):
                    # Extract species ID from beginning of header
                    species_id = line.split('_')[0].replace('>', '')
                    taxid = taxid_map.get(species_id, '0')  # Fallback to 0 if not found
                    # Insert taxID after the first field
                    header_parts = line.strip().split('|')
                    if len(header_parts) >= 3:
                        new_header = f"{header_parts[0].strip()} | taxID:{taxid} | {header_parts[1].strip()} | {header_parts[2].strip()}\n"
                    else:
                        # If malformed, just inject taxID after first pipe
                        new_header = line.strip() + f" | taxID:{taxid}\n"
                    tmp_out.write(new_header)
                else:
                    tmp_out.write(line)

    # Replace original FASTA file
    os.replace(tmp_file, fasta_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <mapping.tsv> <input.fasta>")
        sys.exit(1)

    mapping_file = sys.argv[1]
    fasta_file = sys.argv[2]
    tmp_file = sys.argv[3]

    taxid_map = parse_taxid_mapping(mapping_file)
    update_fasta_headers(fasta_file, tmp_file, taxid_map)