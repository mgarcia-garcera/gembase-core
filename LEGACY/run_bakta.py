import os
import subprocess
import argparse
import sys
import shutil
from GembaseFunctions import run_bakta


def main():
    parser = argparse.ArgumentParser(description="Run Bakta annotation on a FASTA file.")
    
    parser.add_argument("-in", "--input", required=True, help="Input FASTA file")
    parser.add_argument("--db", required=False, help="Path to Bakta database (if not provided, uses $BAKTA_DB)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-f","--force", action="store_true", help="Force overwrite of output directory if it exists")
    
    args = parser.parse_args()
    
    run_bakta(
        fasta_file=args.input,
        output_dir=args.outdir,
        prefix=args.prefix,
        threads=args.threads,
        db_path=args.db,
        force=args.force
    )

if __name__ == "__main__":
    main()