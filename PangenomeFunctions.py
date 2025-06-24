import subprocess
import os

def run_vsearch_clustering(input_fasta, identity=0.5, min_seq_length=50, vsearch_path='vsearch'):
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
    sorted_fasta = f"{file_id}.sort.fasta"
    centroids_fasta = f"{file_id}.centroids.fasta"
    clusters_file = f"{file_id}.clusters.uc"

    try:
        # 1. Sort sequences by length
        subprocess.run([
            vsearch_path,
            '--sortbylength', input_fasta,
            '--fastaout', sorted_fasta,
            '--minseqlength', str(min_seq_length)
        ], check=True)

        # 2. Cluster sequences using small memory algorithm
        subprocess.run([
            vsearch_path,
            '--cluster_smallmem', sorted_fasta,
            '--id', str(identity),
            '--centroids', centroids_fasta,
            '--uc', clusters_file
        ], check=True)

        return {
            'sorted_fasta': sorted_fasta,
            'centroids_fasta': centroids_fasta,
            'clusters_file': clusters_file
        }

    except subprocess.CalledProcessError as e:
        print(f"Error during vsearch execution: {e}")
        return None
