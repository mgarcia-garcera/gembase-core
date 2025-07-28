import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import argparse
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import pandas as pd
import sys
import os

def print_help():
    print("""
Available Clustering Methods:
----------------------------
- ward     : Minimizes the total within-cluster variance (best for Euclidean distances only).
- average  : Uses the average of distances between all pairs in two clusters (UPGMA).
- complete : Uses the maximum distance between elements of each cluster (max linkage).
- single   : Uses the minimum distance between elements of each cluster (min linkage).

Recommended for ANI data: 'average'
""")

def detect_format(input_file):
    with open(input_file) as f:
        header = f.readline().strip().split('\t')
        return 'ANI' in header and 'Ref_file' in header
def load_matrix_format(input_file):
    with open(filename) as f:
        lines = f.readlines()[1:]  # Skip first line
    labels = []
    data = []
    for line in lines:
        parts = line.strip().split('\t')
        labels.append(parts[0])
        data.append([float(x) for x in parts[1:]])
    df = pd.DataFrame(data, index=labels, columns=labels)
    return df

def load_pairwise_format(filename):
    df = pd.read_csv(filename, sep='\t')
    # Use genome names
    ref_names = df['Ref_name'].astype(str)
    query_names = df['Query_name'].astype(str)
    ani = df['ANI'].astype(float)

    # Construct pivot table
    full_df = pd.DataFrame({
        'ref': ref_names,
        'query': query_names,
        'ani': ani
    })
    matrix = full_df.pivot(index='ref', columns='query', values='ani')
    # Fill diagonals and make symmetric
    matrix = matrix.combine_first(matrix.T)
    np.fill_diagonal(matrix.values, 100.0)
    return matrix

def read_ani_matrix(filename):
    with open(filename, "r") as f:
        lines = f.readlines()[1:]  # Skip the first line

    labels = []
    data = []

    for line in lines:
        parts = line.strip().split('\t')
        labels.append(parts[0])
        data.append([float(x) for x in parts[1:]])

    df = pd.DataFrame(data, index=labels, columns=labels)
    return df


def read_distance_matrix(filename):
    labels = []
    data = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            labels.append(parts[0])
            data.append([float(x) for x in parts[1:]])

    matrix = np.array(data)
    return labels, matrix

def validate_distance_matrix(matrix):
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square.")
    if not np.allclose(matrix, matrix.T):
        raise ValueError("Matrix must be symmetric.")
    if not np.all(np.diag(matrix) == 0):
        raise ValueError("Diagonal must be all zeros.")

def plot_dendrogram(labels, matrix, method='ward'):
    condensed = squareform(matrix)
    Z = linkage(condensed, method=method)
    plt.figure(figsize=(8, 4))
    dendrogram(Z, labels=labels)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.show()

def find_centroid(cluster_members, dist_df):
    # Get submatrix for this cluster
    sub = dist_df.loc[cluster_members, cluster_members]
    # Calculate mean distance per genome
    avg_dist = sub.mean(axis=1)
    # Return the genome with minimum average distance
    return avg_dist.idxmin()

def main():
    print("Hierarchical clustering from ANI input\n")
    parser = argparse.ArgumentParser(
        description="""
        Concatenates Fasta files for clustering
        """
    )
    #Parses the fasta headers and provides a GEMBASES identifier, according to AABB.XXX.YYYY
    parser.add_argument("-i", "--input", required=True, help="Input NxN matrix file")
    parser.add_argument("-o", "--output", required=True, help="Output Plot file for plot")
    parser.add_argument("-c", "--centroids", required=True, help="Output Centroid file")
    parser.add_argument("-d", "--distance", default=0.1, help="Distance cutoff to dividing the clusters[0.05-1] . (Default = 0.1)")
    parser.add_argument("-m","--method",required=False, help="Method for the hierarchical clustering")
    parser.add_argument("-f","--format",default="PDF", help="Output Format")
    args = parser.parse_args()

    mthd=args.method
    fmt=args.format
    plot_path = args.output
    centroid_path = args.centroids
    cutoff=args.distance

    if (not args.method) or args.method.lower() == "help":
        print_help()
        return

    if mthd not in ["ward","average","complete","single"]:
        print(f"{args.method} not valid. Use -h for help")
        return 
    
    if fmt.lower() not in ["png", "pdf"]:
        print("Invalid format. Chosen 'pdf' as default.")
        fmt=str("pdf")

    if not os.path.exists(args.input):
        print(f"Input file {args.input} could not be found")
        return
    try:
        is_pairwise = detect_format(args.input)
        if is_pairwise:
            print("Detected long pairwise format")
            df = load_pairwise_format(args.input)
        else:
            print("detected square matrix format")
            df = load_matrix_format(args.input)

        print(f"Loaded matrix of size {df.shape[0]}x{df.shape[1]}")
        dist_df = 100 - df
        dist_df = (dist_df + dist_df.T) / 2  # Ensure symmetry
        np.fill_diagonal(dist_df.values, 0)
        condensed = squareform(dist_df)

        Z = linkage(condensed, method=mthd)

        # Assign clusters based on cut_distance
        cluster_assignments = fcluster(Z, t=cutoff, criterion='distance')
        cluster_dict = {}
        for label, cluster_id in zip(df.index, cluster_assignments):
            cluster_dict.setdefault(cluster_id, []).append(label)
        
        # Output clusters and centroids
        with open(centroid_path, 'w') as f:
            f.write(f"Clusters at distance threshold: {cutoff}\n")
            f.write(f"Total clusters: {len(cluster_dict)}\n\n")
            for cid, members in sorted(cluster_dict.items()):
                centroid = find_centroid(members, dist_df)
                f.write(f"Cluster {cid} ({len(members)} genomes):\n")
                f.write("  Members:\n")
                for m in members:
                    f.write(f"    - {m}\n")
                f.write(f"  Centroid: {centroid}\n\n")

        print(f"Cluster info saved to {centroid_path}")

        plt.figure(figsize=(12, 6))
        dendrogram(Z, labels=df.index.tolist(), leaf_rotation=90)
        plt.title(f"Hierarchical Clustering Dendrogram ({mthd.capitalize()} Linkage)")
        plt.xlabel("Genomes")
        plt.ylabel("Distance (100 - ANI)")
        plt.tight_layout()

        plt.savefig(plot_path)
        print(f"Dendrogram saved to {plot_path}")

        
        
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

