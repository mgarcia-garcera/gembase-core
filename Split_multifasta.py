import csv
from Bio import SeqIO
from collections import defaultdict
import os
import gzip

def is_gzipped(filepath):
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def open_maybe_gzipped(filepath):
    if is_gzipped(filepath):
        return gzip.open(filepath, 'rt')  # text mode
    else:
        return open(filepath, 'r')

def load_metadata(tsv_file):
    accession_to_info = {}
    with open(tsv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            acc = row['AccessionID'].strip()
            gembase_id = row['GembasesID'].strip()
            taxid = row['TaxID'].strip()
            accession_to_info[acc] = {
                'GembasesID': gembase_id,
                'TaxID': taxid
            }
    return accession_to_info

def split_fasta_by_gembaseid(fasta_file, metadata_tsv,output_dir):
    accession_map = load_metadata(metadata_tsv)

    output_handles = {}  # GembasesID -> file handle

    with open_maybe_gzipped(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header_parts = record.description.split()
            accession = header_parts[0]

            if accession not in accession_map:
                print(f"Warning: Accession {accession} not found in metadata. Skipping.")
                continue

            info = accession_map[accession]
            gembase_id = info['GembasesID']
            taxid = info['TaxID']


        # Modify header
            new_header = f"{accession}|taxID:{taxid}|{' '.join(header_parts[1:])}"
            record.description = ""
            record.id = new_header
            record.name = new_header
        # Create output file if not open yet
            output_file = os.path.join(output_dir,f"{gembase_id}.fna")
            with open(output_file, "a") as out_handle:
                SeqIO.write(record,out_handle, "fasta")


    print("Done. Sequences written into individual GembasesID.fna files.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Split a multifasta file into GembasesID.fna files based on metadata TSV.")
    parser.add_argument("-i", help="Input multi-FASTA file")
    parser.add_argument("-m", help="Metadata TSV file (AccessionID<TAB>GembasesID<TAB>TaxID)")
    parser.add_argument("-o", help="Output directory")
    args = parser.parse_args()

    split_fasta_by_gembaseid(args.i, args.m,args.o)