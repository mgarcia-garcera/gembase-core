import subprocess
import os
from pathlib import Path
import tempfile
from Bio import SeqIO
import time
from functools import wraps
import shutil
import sys
from collections import defaultdict
import gzip
def timed(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        duration = end - start
        print(f"⏱️ Function '{func.__name__}' took {duration:.2f} seconds")
        return result
    return wrapper

def extract_fam_file(fam_file, final_fam):
    shutil.copy(fam_file,final_fam)
    print(f"Copied FAM file to: {final_fam}")
    return 

def concatenate_species_fastas(identifier, input_folder, tmpdirname):
    """
    Concatenates all FASTA files in input_folder whose names contain the species identifier into a single file.
    
    Parameters:
        identifier (str): Species identifier (e.g. "AABB.XXX")
        input_folder (str): Path to the folder containing input FASTA files
        output_file (str): Name/path of the output concatenated FASTA file (default: 'tmp.fasta')
    
    Returns:
        int: Number of files concatenated
    """
    input_path = Path(input_folder)
    matching_files = [f for f in input_path.glob("*") if identifier in f.name]

    if not matching_files:
        print(f"❌ No files found with identifier '{identifier}' in {input_folder}")
        return 0

    output_file = os.path.join(tmpdirname, f"{identifier}.tmp.fna")
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for fasta_file in matching_files:
            with open(fasta_file, 'r', encoding='utf-8') as infile:
                outfile.write(infile.read())
                if not infile.read().endswith('\n'):
                    outfile.write('\n')  # Ensure newline between files

    print(f"✅ Concatenated {len(matching_files)} files into {output_file}")
    return output_file

def run_cmd(cmd, log_path):
    with open(log_path, "a") as log_file:
        subprocess.run(
            cmd,
            stdout=log_file,
            stderr=log_file,
            check=True
        )
def rename_fasta_headers(input_fasta, output_fasta, species_id, taxids):
    """
    Renames FASTA headers to:
        >{species_id}_FAMXXXXXX0 | description | strain_id

    Args:
        input_fasta (str or Path): Input FASTA file with centroid sequences.
        output_fasta (str or Path): Output FASTA file with renamed headers.
        species_id (str): Species prefix (e.g. "BAAL.001").
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    renamed_records = []
    counter = 1

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Step 1: Get full original ID
        original_id = record.id  # e.g., BAAL.001.0028.d01_000001

        # Step 2: Extract strain ID (first 4 dot-separated fields)
        strain_parts = original_id.split(".")
        strain_id = ".".join(strain_parts[:4]) if len(strain_parts) >= 4 else original_id.split("_")[0]

        # Step 3: Extract description after coordinate fields (starts at field 5+)
        desc_fields = record.description.split()
        description = " ".join(desc_fields[5:]) if len(desc_fields) > 5 else ""

        # Step 4: Build new ID and description
        fam_id = f"{species_id}_FAM{counter:07d}0"
        record.id = fam_id
        if taxids:
            tax_id = taxids[species_id]
            record.description = f"{fam_id} | taxID:{tax_id} | {description} | {strain_id}"
        else:
            record.description = f"{fam_id} | {description} | {strain_id}"
                
        renamed_records.append(record)
        counter += 1

    cnt = len(renamed_records)
    # Write new FASTA
    with open(output_fasta, "w") as out_f:
        SeqIO.write(renamed_records, out_f, "fasta")

    print(f"✅ Renamed {cnt} records → {output_fasta}")
    return cnt


def count_files_with_speciesId(folder_path, keyword):
    folder = Path(folder_path)
    matching_files = [f for f in folder.iterdir() if f.is_file() and keyword in f.name]
    return len(matching_files)

def get_value_by_key(filename, keyword, delimiter="\t"):
    """
    Reads a 2-column file and returns the value in column 2 that matches the keyword in column 1.

    Args:
        filename (str): Path to the input file.
        keyword (str): Value to search for in column 1.
        delimiter (str or None): Field separator (default: any whitespace).

    Returns:
        str or None: Value from column 2, or None if not found.
    """
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip() == "":
                continue  # skip empty lines
            parts = line.strip().split(delimiter)
            if len(parts) < 2:
                continue  # skip malformed lines
            if parts[0] == keyword:
                return parts[1]
    return None  # Not found 

def run_vsearch_clustering(input_fasta, output_folder, identity=0.8, min_seq_length=50, threads=24, min_cov=0.9, vsearch_path='vsearch', log="output.log"):
    """
    Runs vsearch to sort sequences by length and perform clustering.
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        identity (float): Clustering identity threshold (default=0.5).
        min_seq_length (int): Minimum sequence length to retain (default=50).
        vsearch_path (str): Path to the vsearch executable (default='vsearch' assumes it's in PATH).
    
    Returns:
        dict: Paths to the generated output files.
    """
    file_id = os.path.splitext(os.path.basename(input_fasta))[0]
    sorted_fasta = os.path.join(output_folder, f"{file_id}.sort.fasta")
    centroids_fasta = os.path.join(output_folder, f"{file_id}.centroids.fasta")
    clusters_file = os.path.join(output_folder, f"{file_id}.clusters.uc")

    try:
        # 1. Sort sequences by length
        cmd1 = [
            vsearch_path,
            '--sortbylength', input_fasta,
            '--output', sorted_fasta,
            '--minseqlength', str(min_seq_length)
        ]
        run_cmd(cmd1, log)

        # 2. Cluster sequences using small memory algorithm
        cmd2 = [
            vsearch_path,
            '--cluster_smallmem', sorted_fasta,
            '--id', str(identity),
            '--centroids', centroids_fasta,
            '--uc', clusters_file,
            '--minqt', str(min_cov),
            '--threads', str(threads)
        ]
        run_cmd(cmd2, log)

        return centroids_fasta, clusters_file

    except subprocess.CalledProcessError as e:
        print(f"Error during vsearch execution: {e}")
        return None

def run_mmseq2_clustering(input_fasta, output_dir,
                          min_identity=0.9, min_coverage=0.8, force=False, log="output.log"):
    """
    Perform MMseqs2 clustering on nucleotide sequences and return a mapping table.

    Args:
        input_fasta (str or Path): Input concatenated nucleotide FASTA file.
        output_dir (str or Path): Directory for intermediate and output files.
        min_identity (float): Minimum sequence identity (0–1).
        min_coverage (float): Minimum coverage threshold (0–1).
        force (bool): Re-run even if output exists.

    Returns:
        Path: Path to TSV file with columns [representative, member]
    """
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_db = output_dir / "clusterDB"
    cluster_result = output_dir / "clusterRes"
    mapping_tsv = output_dir / "cluster_mapping.tsv"
    tmp_dir = tempfile.mkdtemp()

    if mapping_tsv.exists() and not force:
        print(f"⚠️ {mapping_tsv} exists. Skipping clustering (use force=True to override).")
        return mapping_tsv

    # Step 1: Create DB
    cmd1 = ["mmseqs", "createdb", str(input_fasta), str(cluster_db)]
    run_cmd(cmd1, log)

    # Step 2: Cluster
    cmd2 = [
        "mmseqs", "cluster",
        str(cluster_db), str(cluster_result), tmp_dir,
        "--min-seq-id", str(min_identity),
        "-c", str(min_coverage),
        "--cov-mode", "0",
        "--cluster-mode", "2",
        "--seq-id-mode", "1"
    ]
    run_cmd(cmd2, log)

    # Step 3: Extract mapping
    cmd3 = [
        "mmseqs", "createtsv",
        str(cluster_db), str(cluster_db), str(cluster_result),
        str(mapping_tsv)
    ]
    run_cmd(cmd3, log)
    print(f"✅ MMseqs2 clustering complete.")
    print(f"  Cluster mapping TSV: {mapping_tsv}")
    return mapping_tsv

def extract_representative_fasta(tsv_path, original_fasta, output_fasta):
    """
    Extracts representative sequences based on MMseqs2 cluster TSV.

    Args:
        tsv_path (str or Path): Path to TSV with [rep_id, member_id] per line.
        original_fasta (str or Path): Original nucleotide FASTA file.
        output_fasta (str or Path): Where to write representative sequences.
    """
    tsv_path = Path(tsv_path)
    original_fasta = Path(original_fasta)
    output_fasta = Path(output_fasta)

    # Step 1: Parse TSV and collect all representative IDs
    rep_ids = set()
    with open(tsv_path, "r") as f:
        for line in f:
            rep_id, _ = line.strip().split()
            rep_ids.add(rep_id)

    # Step 2: Read original FASTA and extract matching IDs
    with open(original_fasta, "r") as infile, open(output_fasta, "w") as out:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in rep_ids:
                SeqIO.write(record, out, "fasta")

    print(f"✅ Extracted {len(rep_ids)} representative sequences → {output_fasta}")
    return output_fasta

def run_cd_hit_est_clustering(input_fasta, output_dir, identity=0.95, coverage=0.8, memory=0, word_size=None, threads=4, force=False, log="output.log"):
    """
    Cluster nucleotide sequences using CD-HIT-EST.

    Args:
        input_fasta (str or Path): Path to input FASTA file (nucleotide).
        output_dir (str or Path): Directory to save outputs.
        identity (float): Sequence identity threshold (0.0 - 1.0, e.g., 0.95).
        coverage (float): Minimum coverage for shorter sequence (0.0 - 1.0).
        word_size (int): Optional; word size (default is inferred by cd-hit-est).
        threads (int): Number of threads to use.
        force (bool): If False, skip if output exists.

    Outputs:
        - output_dir/clustered.fasta
        - output_dir/clustered.fasta.clstr
    """
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = output_dir / "clustered.fasta"
    clstr_file = str(output_fasta) + ".clstr"
    if output_fasta.exists() and not force:
        print(f"⚠️ Output {output_fasta} exists. Skipping CD-HIT-EST (use force=True to override).")
        return output_fasta, clstr_file
    cmd = [
        "cd-hit-est",
        "-i", str(input_fasta),
        "-o", str(output_fasta),
        "-c", str(identity),
        "-aS", str(coverage),
        "-T", str(threads), 
        "-M", str(memory)
    ]

    if word_size:
        cmd.extend(["-n", str(word_size)])

    print(f" Running CD-HIT-EST with identity={identity}, coverage={coverage}")
    run_cmd(cmd, log)
    print(f"✅ CD-HIT-EST clustering complete.")
    print(f"  - Representatives: {output_fasta}")
    print(f"  - Cluster info:    {clstr_file}")

    return output_fasta, clstr_file

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


class UnionFind:
    def __init__(self):
        self.parent = {}
    
    def find(self, item):
        if item not in self.parent:
            self.parent[item] = item
        if self.parent[item] !=item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]
    
    def union(self, a, b):
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a != root_b:
            self.parent[root_b] = root_a


def process_and_cluster(input_path, output_path,
                        identity_threshold=95.0,
                        coverage1_threshold=70.0,
                        coverage2_threshold=70.0):
    uf = UnionFind()
    edges = defaultdict(list)

    print("First pass: Reading and building clusters...")

    with open(input_path, 'r') as infile:
        for line_num, line in enumerate(infile, 1):
            try:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 7:
                    continue

                identity = float(fields[2])
                coverage1 = float(fields[3])
                coverage2 = float(fields[4])

                if (identity < identity_threshold or
                    coverage1 < coverage1_threshold or
                    coverage2 < coverage2_threshold):
                    continue  # skip low-quality comparisons

                acc1 = fields[5].split()[0]
                acc2 = fields[6].split()[0]

                # Build connections
                uf.union(acc1, acc2)
                edges[acc1].append((acc2, identity))
                edges[acc2].append((acc1, identity))

                if line_num % 1_000_000 == 0:
                    print(f"Processed {line_num} lines...")

            except Exception as e:
                print(f"Error at line {line_num}: {e}", file=sys.stderr)
                continue

    print("Second pass: Assembling clusters...")

    clusters = defaultdict(set)
    for accession in uf.parent:
        root = uf.find(accession)
        clusters[root].add(accession)

    print(f"Found {len(clusters)} clusters.")

    with open(output_path, 'w') as outfile:
        for idx, (cluster_root, members) in enumerate(clusters.items(), 1):
            clstr_id = f"CLSTR_{idx:04d}"
            entries = []
            seen = set()

            for acc in members:
                for neighbor, ident in edges[acc]:
                    if neighbor in members and (acc, neighbor) not in seen and (neighbor, acc) not in seen:
                        entries.append(f"{acc}:{ident:.2f}")
                        entries.append(f"{neighbor}:{ident:.2f}")
                        seen.add((acc, neighbor))

            # Remove duplicates and preserve order
            unique_entries = list(dict.fromkeys(entries))
            outfile.write(f"{clstr_id}\t" + "\t".join(unique_entries) + "\n")


def extract_accessions_from_gbff(gbff_path, output_path=None):
    accessions = []

    with gzip.open(gbff_path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            acc = record.id.split('.')[0]  # Remove version if needed
            accessions.append(acc)

    return accessions

def extract_centroid_sequences(clustering_file, fasta_file, output_fasta):
    """
    Extracts sequences from a FASTA file whose headers match the centroids
    (unique values from column 1) in the clustering file.
    
    Parameters:
        clustering_file (str): Path to the clustering file.
        fasta_file (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file containing centroid sequences.
    """
    # Extract unique centroids from column 1
    centroids = set()
    with open(clustering_file, "r") as f:
        for line in f:
            if line.strip():
                centroid = line.strip().split()[0]
                centroids.add(centroid)

    # Parse FASTA and extract matching records
    records = SeqIO.parse(fasta_file, "fasta")
    matching_records = (record for record in records if record.id in centroids)

    # Write the results
    count = SeqIO.write(matching_records, output_fasta, "fasta")
    print(f"Extracted {count} centroid sequences to '{output_fasta}'")

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


