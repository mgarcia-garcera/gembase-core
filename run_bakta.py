import argparse
import subprocess
import os
import sys

def run_bakta(fasta_path, output_dir, path2db, threads=4, prefix=None):
    if not os.path.isfile(fasta_path):
        print(f"[ERROR] FASTA file not found: {fasta_path}")
        sys.exit(1)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    cmd = [
        "bakta",
        "--db", "/path/to/bakta/db",  # UPDATE this path as appropriate
        "--output", output_dir,
        "--prefix", prefix if prefix else os.path.splitext(os.path.basename(fasta_path))[0],
        "--force",
        "--threads", str(threads),
        fasta_path
    ]

    print(f"[INFO] Running Bakta:\n{' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] Bakta finished. GenBank file in: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Bakta failed: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run Bakta on a FASTA file and generate a GenBank file.")
    parser.add_argument("fasta", help="Path to the .renamed.fasta file")
    parser.add_argument("-o", "--output", default="bakta_output", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for Bakta")
    parser.add_argument("--prefix", help="Prefix for Bakta output files")

    args = parser.parse_args()
    run_bakta(args.fasta, args.output, args.threads, args.prefix)

if __name__ == "__main__":
    main()
