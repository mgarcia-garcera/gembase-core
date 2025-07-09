import argparse
import os
import tempfile
from pathlib import Path
import subprocess
from Bio import SeqIO


def concatenate_species_fastas(identifier, input_folder, output_file="tmp.fasta"):
    input_folder = Path(input_folder)
    output_path = Path(output_file)
    with open(output_path, "w") as out_handle:
        for fasta_file in input_folder.glob("*.fasta"):
            if identifier in fasta_file.name:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    SeqIO.write(record, out_handle, "fasta")
    print(f"✅ Concatenated FASTAs written to {output_file}")


def run_mmseq2_clustering(identifier, input_fasta, output_dir, min_identity=0.8, min_coverage=0.8, force=False):
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_db = output_dir / "clusterDB"
    tmp_dir = tempfile.mkdtemp()

    # Step 1. Create MMseqs2 database
    subprocess.run(["mmseqs", "createdb", str(input_fasta), str(cluster_db)], check=True)

    # Step 2. Cluster sequences
    cluster_result = output_dir / "clusterRes"
    subprocess.run([
        "mmseqs", "cluster",
        str(cluster_db), str(cluster_result), tmp_dir,
        "--min-seq-id", str(min_identity),
        "-c", str(min_coverage),
        "--cov-mode", "0"
    ], check=True)

    # Step 3. Extract representative sequences
    rep_file_name = f"{identifier}.reps.fasta"
    rep_fasta = output_dir / rep_file_name

    if not rep_fasta.exists() or force:
        subprocess.run(["mmseqs", "createseqfiledb",
                        str(cluster_db), str(cluster_result), str(output_dir / "repDB")], check=True)
        subprocess.run(["mmseqs", "result2flat",
                        str(cluster_db), str(cluster_db), str(output_dir / "repDB"), str(rep_fasta)], check=True)
    else:
        print(f"⚠️ Representative file {rep_fasta} already exists. Use --force to overwrite.")

    # Step 4. Full cluster mappings
    tsv_path = output_dir / "cluster_mapping.tsv"
    subprocess.run([
        "mmseqs", "createtsv",
        str(cluster_db), str(cluster_db), str(cluster_result), str(tsv_path)
    ], check=True)

    # Step 5a. Generate FAM names
    fam_file = output_dir / "clusters_with_fam_names.tsv"
    fam_map = {}
    cluster_counter = 1

    with open(tsv_path, "r") as tsv_in, open(fam_file, "w") as fam_out:
        for line in tsv_in:
            rep, member = line.strip().split()
            if rep not in fam_map:
                fam_name = f"FAM{cluster_counter:07d}"
                fam_map[rep] = fam_name
                cluster_counter += 1
            fam_out.write(f"{fam_map[rep]}\t{member}\n")

    # Step 5b. Alignment statistics
    aln_tsv = output_dir / "cluster_alignments.tsv"
    subprocess.run([
        "mmseqs", "align",
        str(cluster_db), str(cluster_db), str(output_dir / "alignmentDB"), tmp_dir,
        "--min-seq-id", str(min_identity),
        "-c", str(min_coverage),
        "--cov-mode", "0"
    ], check=True)
    subprocess.run([
        "mmseqs", "createtsv",
        str(cluster_db), str(cluster_db), str(output_dir / "alignmentDB"), str(aln_tsv)
    ], check=True)

    print(f"✅ MMseqs2 clustering completed.")
    print(f"Representatives: {rep_fasta}")
    print(f"Cluster members: {fam_file}")
    print(f"Identity/Coverage: {aln_tsv}")

    return rep_fasta, fam_file


def rename_sequences_with_fam(rep_fasta, fam_file, identifier, output_file):
    fam_dict = {}
    for line in open(fam_file):
        fam, member = line.strip().split("\t")
        fam_dict[member] = fam

    with open(rep_fasta) as in_fh, open(output_file, "w") as out_fh:
        for record in SeqIO.parse(in_fh, "fasta"):
            orig_id = record.id
            fam = fam_dict.get(orig_id, "FAM0000000")
            record.id = f"{identifier}.{fam}0"
            record.description = f"{record.description} |ORIG: {orig_id}|"
            SeqIO.write(record, out_fh, "fasta")
    print(f"✅ Sequences renamed and saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Cluster species FASTAs with MMseqs2")
    parser.add_argument("-s", "--species", required=True, help="Species identifier (e.g., AABB.XXX)")
    parser.add_argument("-i", "--input_folder", required=True, help="Path to folder with FASTA files")
    parser.add_argument("-m", "--method", required=False, default="mmseqs", help="Clustering method (only 'mmseqs' supported)")
    parser.add_argument("--force", action="store_true", help="Force rerun if outputs exist")
    parser.add_argument("-o", "--output_folder", default="mmseqs_output", help="Output directory")
    args = parser.parse_args()

    if args.method != "mmseqs":
        raise ValueError("Only 'mmseqs' method is currently supported.")

    tmp_fasta = "tmp.fasta"
    concatenate_species_fastas(args.species, args.input_folder, tmp_fasta)
    rep_fasta, fam_file = run_mmseq2_clustering(
        args.species, tmp_fasta, args.output_folder, force=args.force
    )
    final_output = Path(args.output_folder) / f"{args.species}.renamed_reps.fasta"
    rename_sequences_with_fam(rep_fasta, fam_file, args.species, final_output)


if __name__ == "__main__":
    main()