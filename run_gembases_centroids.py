import argparse
import os
import re
import subprocess


def load_accession_map(accession_file):
    acc_map = {}
    with open(accession_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                acc, path = parts
                acc_map[acc] = path 
    return acc_map


def parse_centroid_accessions(centroids_file):
    centroids = []
    with open(centroids_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Centroid:"):
                
                match = re.search(r'Centroid:\s([\w\.\d]+)', line)  
                if match: 
                    centroids.append(match.group(1))
                    
    return centroids

def run_builder_script(gbff_file, reference, history, output, email,log):
    cmd = [
        "python", "/scratch/hdd3/mgg/notebooks/labbooks/Marc/Gembases2NextflowGembases/gembase-core/gembases_builderGBFF.py", 
        "-i", gbff_file,
        "-r", reference,
        "-hi", history,
        "-o", output,
        "-e", email,
        "-l", log
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Processed {gbff_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {gbff_file}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Run gembases_builderGBFF.py on centroid genomes.")
    parser.add_argument("-a", "--accessions", required=True, help="TSV file from get_accessions.py")
    parser.add_argument("-c", "--centroids", required=True, help="Centroids file with clusters")
    parser.add_argument("-r", "--reference", required=True, help="Reference file path")
    parser.add_argument("-hi", "--history", required=True, help="History file path")
    parser.add_argument("-o", "--output", required=True, help="Output folder path")
    parser.add_argument("-e", "--email", required=True, help="Email address")
    parser.add_argument("-l","--log",help="logfolder")
    args = parser.parse_args()

    accession_map = load_accession_map(args.accessions)
    centroid_accessions = parse_centroid_accessions(args.centroids)

    try:
        os.path.exists(args.accessions)
    except Exception as e:
        f"❌ ERROR: {e}"
    
    try:
        os.path.exists(args.centroids)
    except Exception as e:
        f"❌ ERROR: {e}"
    if os.path.exists(args.accessions) and os.path.exists(args.centroids):
        for accession in centroid_accessions:
            gbff_file = accession_map.get(accession)
            logfile = os.path.join(args.log, f"{accession}.log")
            if gbff_file:
                run_builder_script(gbff_file, args.reference, args.history, args.output, args.email, logfile)
            else:
                print(f"⚠️ Accession not found in map: {accession}")

if __name__ == "__main__":
    main()