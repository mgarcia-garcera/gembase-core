import os
import gzip
from Bio import Entrez, SeqIO


Entrez.email = "your.email@example.com"  # <-- Replace this!

def is_gzipped(filepath):
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def open_maybe_gzipped(filepath):
    return gzip.open(filepath, 'rt') if is_gzipped(filepath) else open(filepath, 'r')

def fetch_taxid(accession):
    try:
        handle = Entrez.esummary(db="nucleotide", id=accession, retmode="xml")
        summary = Entrez.read(handle)
        handle.close()

        docsum = summary[0]

        # Correctly unwrap IntegerElement and convert to string
        taxid_raw = docsum.get("TaxId", None)
        taxid = str(int(taxid_raw)) if taxid_raw is not None else "NA"
        return taxid

    except Exception as e:
        print(f"Failed to fetch taxid for {accession}: {e}")
        return "NA"
    
def update_fasta_headers(input_fasta, output_fasta):
    seen = {}
    with open_maybe_gzipped(input_fasta) as handle, open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(handle, "fasta"):
            parts = record.description.split()
            accession = parts[0]

            # Avoid repeated API calls
            if accession not in seen:
                seen[accession] = fetch_taxid(accession)

            taxid = seen[accession]

            # Modify header
            new_header = f"{accession}|taxID:{taxid}|{' '.join(parts[1:])}"
            record.id = new_header
            record.name = new_header
            record.description = ""

            SeqIO.write(record, out_handle, "fasta")

    print(f"Finished writing updated FASTA to {output_fasta}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Add NCBI TaxID to FASTA headers using Entrez.")
    parser.add_argument("-i", help="Input multi-FASTA file (can be gzipped)")
    parser.add_argument("-o", help="Output FASTA file with |taxID:xxxxxx| in headers")
    args = parser.parse_args()

    update_fasta_headers(args.i, args.o)