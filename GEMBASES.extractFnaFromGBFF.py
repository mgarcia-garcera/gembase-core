import argparse
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
import os
from multiprocessing import Pool, cpu_count


def extract_fasta_from_gbff(args):
    gbff_path, output_folder = args

    base_name = os.path.basename(gbff_path).replace(".gbff.gz",".fna.gz")
    output_path = os.path.join(output_folder, base_name)

    try: 
        with gzip.open(gbff_path, "rt") as gbff_handle:
            records = list(SeqIO.parse(gbff_handle,"genbank"))
    except Exception as e: 
        print(f"❌ Error reading {gbff_path}:{e}")
        return
    
    if not records:
        print(f"❌ No records found in {gbff_path}")

    try:
        with gzip.open(output_path, "wt") as fna_handle:
            fasta_writer = FastaIO.FastaWriter(fna_handle, wrap=60)
            for record in records:
                fasta_writer.write_record(SeqRecord(record.seq, id=record.id, description=record.description))
    except Exception as e:
        print(f"❌ Error writing {output_path}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Extract nucleotide FASTA from a .gbff.gz file")
    parser.add_argument("-i", "--input", required=True, help="Path to .gbff.gz files")
    parser.add_argument("-o","--output", help="Optional output directory to print .fna.gz")
    parser.add_argument("-p", "--processes", type=int, default=10, help="Number of parallel processes")
    args = parser.parse_args()

    input_folder = os.path.abspath(args.input)
    output_folder = os.path.abspath(args.output) if args.output else input_folder

    num_processes = args.processes

    gbff_files =[os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".gbff.gz")]

    if not gbff_files:
        print(f"❌ No .gbff.gz files found in input folder {input_folder}")
        return

    print(f"  Found {len(gbff_files)} files. Starting with {num_processes} process(es)...")

    with Pool(processes=num_processes) as pool:
        pool.map(extract_fasta_from_gbff,[(f, output_folder) for f in gbff_files])

if __name__ == "__main__":
    main()