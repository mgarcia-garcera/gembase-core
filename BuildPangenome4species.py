import argparse
import os
import tempfile
from pathlib import Path
import subprocess
from Bio import SeqIO
import PangenomeFunctions as PF
import GembaseFunctionsGBFF as GF
import tempfile

@PF.timed 
def main():
    parser = argparse.ArgumentParser(description="Cluster species FASTAs with MMseqs2")
    parser.add_argument("-s", "--species", required=True, help="Species identifier (e.g., AABB.XXX)")
    parser.add_argument("-f", "--input_folder", required=True, help="Path to folder with FASTA files")
    parser.add_argument("-m", "--method", required=False, default="mmseqs", help="Clustering method")
    parser.add_argument("--force", action="store_true", help="Force rerun if outputs exist")
    parser.add_argument("-o", "--output_folder", default="mmseqs_output", help="Output directory")
    parser.add_argument("-i","--min_identity", required=False, default=0.8, help="Identity threshold")
    parser.add_argument("-c","--min_coverage", required=False, default=0.8, help="Coverage threshold")
    parser.add_argument("-t","--threads", required=False, default=24, help="Number of threads")
    parser.add_argument("-l", "--log", required=False, help="Clustering Log output")
    args = parser.parse_args()

    tmpdir = tempfile.TemporaryDirectory()
    print(f"I have created a temporary folder in {tmpdir.name}")

    if not args.log:
        logfile = os.path.join(args.output_folder,f"{args.species}.{args.method}.log")
        if not os.path.exists(logfile):
            open(logfile, 'a').close()
    else:
        logfile = args.log
        if not os.path.exists(logfile):
            open(logfile, 'a').close()


    tmp_fasta = PF.concatenate_species_fastas(args.species, args.input_folder, tmpdir.name)

    if args.method == 'mmseqs':
        fam_file = PF.run_mmseq2_clustering(tmp_fasta,tmpdir.name, min_identity=args.min_identity, min_coverage=args.min_coverage, force=True, log=logfile)
        rep_fasta = PF.extract_representative_fasta(tsv_path=fam_file, original_fasta=tmp_fasta,output_fasta=f"{tmpdir.name}{args.species}.reps.fna")
    elif args.method == 'cdhit':
        rep_fasta, fam_file = PF.run_cd_hit_est_clustering(tmp_fasta,tmpdir.name,identity=args.min_identity, coverage=args.min_coverage, memory=0, threads=24, force=True, log=logfile)
    elif args.method == 'vsearch':
        rep_fasta, fam_file = PF.run_vsearch_clustering(tmp_fasta, tmpdir.name, identity=args.min_identity, min_seq_length=50, min_cov=args.min_coverage, vsearch_path='vsearch',log=logfile)

    # Modify the fasta to include information about the species:
    final_fasta_name = f"{args.species}.pangenome.{args.method}.i{args.min_identity}.c{args.min_coverage}.fna"
    final_fam_name = f"{args.species}.pangenome.{args.method}.i{args.min_identity}.c{args.min_coverage}.txt"
    final_fasta = os.path.join(args.output_folder, final_fasta_name)
    final_fam = os.path.join(args.output_folder, final_fam_name)
    PF.rename_fasta_headers(rep_fasta, final_fasta, args.species)
    cnt = PF.extract_fam_file(fam_file, final_fam)
#    final_output = Path(args.output_folder) / f"{args.species}.renamed_reps.fasta"
#    PF.rename_sequences_with_fam(rep_fasta, fam_file, args.species, final_output)
    
    tmpdir.cleanup()

    print(f"✅ Pangenome building successfully finished → {final_fasta}")

if __name__ == "__main__":
    main()