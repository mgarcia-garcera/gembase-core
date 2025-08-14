import sys
import os
import tempfile
import argparse

def fix_kraken_headers(fasta_path):
    # Create a temporary output file
    with open(fasta_path, 'r') as infile, open("test2.fna", 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Example: >ABC123.1¦taxID:1206811¦blabla
                parts = line.strip().split('¦')
                if len(parts) >= 2 and parts[1].startswith("taxID:"):
                    seq_id = parts[0][1:]  # remove '>'
                    taxid = parts[1].split(':')[1]
                    description = ' '.join(parts[2:]) if len(parts) > 2 else ''
                    new_header = f">{seq_id}|kraken:taxid|{taxid} {description}".strip()
                    outfile.write(new_header + '\n')
                else:
                    # Write original if it doesn't match the expected format
                    outfile.write(line)
            else:
                outfile.write(line)

    # Replace the original file with the modified one
    os.replace("test2.fna", fasta_path)
    print(f"Updated headers in: {fasta_path}")

# --- Usage ---
# Replace 'your_file.fasta' with your filename or pass via command line
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fix_kraken_headers.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    fix_kraken_headers(fasta_file)