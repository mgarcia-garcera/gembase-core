import argparse
import os
import tempfile
from pathlib import Path
import subprocess
from Bio import SeqIO
import PangenomeFunctions as PF
import GembaseFunctionsGBFF as GF
import tempfile
import shutil
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
    parser.add_argument("--tmp", required=False, help="Temporary folder location: For long clustering, it is needed to provide a local folder. Otherwise the script crashes")
    parser.add_argument("--tax")
    args = parser.parse_args()


    if args.tax:
        taxIDs = PF.parse_taxid_mapping(args.tax)

    if not args.tmp:
        tmpdir = tempfile.TemporaryDirectory()
        tmpdirname = tmpdir.name
        
    else:
        tmpdirname = args.tmp
        if not os.path.exists(tmpdirname):
            os.makedirs(tmpdirname)
        elif os.path.exists(tmpdirname):
            shutil.rmtree(tmpdirname)
            os.makedirs(tmpdirname)

        
    print(f"I have created a temporary folder in {tmpdirname}")

    if not args.log:
        logfile = os.path.join(args.output_folder,f"{args.species}.{args.method}.log")
        if not os.path.exists(logfile):
            open(logfile, 'a').close()
    else:
        logfile = args.log
        if not os.path.exists(logfile):
            open(logfile, 'a').close()


    final_fasta_name = f"{args.species}.pangenome.{args.method}.i{args.min_identity}.c{args.min_coverage}.fna"
    final_fam_name = f"{args.species}.pangenome.{args.method}.i{args.min_identity}.c{args.min_coverage}.txt"
    final_fasta = os.path.join(args.output_folder, final_fasta_name)
    final_fam = os.path.join(args.output_folder, final_fam_name)
    if not os.path.exists(final_fasta):
        tmp_fasta = PF.concatenate_species_fastas(args.species, args.input_folder, tmpdirname)
        if os.path.exists(final_fam):
            PF.extract_centroid_sequences(final_fam, tmp_fasta, output_fasta=f"{tmpdirname}{args.species}.reps.fna")
        else:
            if args.method == 'mmseqs':
                fam_file = PF.run_mmseq2_clustering(tmp_fasta,tmpdirname, min_identity=args.min_identity, min_coverage=args.min_coverage, force=True, log=logfile)
                rep_fasta = PF.extract_representative_fasta(tsv_path=fam_file, original_fasta=tmp_fasta,output_fasta=f"{tmpdirname}{args.species}.reps.fna")
            elif args.method == 'mmseqs2':
                fam_file = PF.run_mmseq2_clustering(tmp_fasta,tmpdirname, min_identity=args.min_identity, min_coverage=args.min_coverage, force=True, log=logfile)
                rep_fasta = PF.extract_representative_fasta(tsv_path=fam_file, original_fasta=tmp_fasta,output_fasta=f"{tmpdirname}{args.species}.reps.fna")
            elif args.method == 'cdhit':
                rep_fasta, fam_file = PF.run_cd_hit_est_clustering(tmp_fasta,tmpdirname,identity=args.min_identity, coverage=args.min_coverage, memory=0, threads=24, force=True, log=logfile)
            elif args.method == 'vsearch':
                rep_fasta, fam_file = PF.run_vsearch_clustering(tmp_fasta, tmpdirname, identity=args.min_identity, min_seq_length=50, min_cov=args.min_coverage, vsearch_path='vsearch',log=logfile)
            else:
                print(f"{args.method} is not known. Running MMseqs2 instead")
                fam_file = PF.run_mmseq2_clustering(tmp_fasta,tmpdirname, min_identity=args.min_identity, min_coverage=args.min_coverage, force=True, log=logfile)
                rep_fasta = PF.extract_representative_fasta(tsv_path=fam_file, original_fasta=tmp_fasta,output_fasta=f"{tmpdirname}{args.species}.reps.fna")

        # Modify the fasta to include information about the species:

        PF.rename_fasta_headers(rep_fasta, final_fasta, args.species, taxIDs)
        cnt = PF.extract_fam_file(fam_file, final_fam)
#    final_output = Path(args.output_folder) / f"{args.species}.renamed_reps.fasta"
#    PF.rename_sequences_with_fam(rep_fasta, fam_file, args.species, final_output)
    else:
        print(f"Pangenome for {args.species} already performed")
    if not args.tmp:    
        tmpdir.cleanup()
    else:
        shutil.rmtree(tmpdirname)

    print(f"✅ Pangenome building successfully finished → {final_fasta}")

if __name__ == "__main__":
    main()