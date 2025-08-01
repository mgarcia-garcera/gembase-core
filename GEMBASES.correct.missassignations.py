import os
import re
import argparse
import shutil
import gzip

from collections import defaultdict

def read_corrections(corrections_file):
    corrections = []
    with open(corrections_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            assigned_id = cols[0].strip()
            expected_id = cols[2].strip()
            location = cols[3].strip()
            corrections.append((assigned_id, expected_id, location))
    return corrections

def find_matching_genome_ids(gembase_root, species_id):
    """Return sorted list of genome IDs like STDS.001.0001, from files in subfolders."""
    pattern = re.compile(rf"({re.escape(species_id)}\.\d{{4}})")
    genome_ids = set()
    for subdir in ['Replicons', 'Genes', 'Proteins', 'RNA', 'LSTINFO']:
        folder = os.path.join(gembase_root, subdir)
        if not os.path.isdir(folder):
            continue
        for root, _, files in os.walk(folder):
            for file in files:
                match = pattern.search(file)
                if match:
                    genome_ids.add(match.group(1))
    return sorted(genome_ids)

def get_last_assigned_number(gembase_root, correct_species_id):
    """Find the highest NNNN in correct_species_id.NNNN across all subfolders."""
    pattern = re.compile(rf"{re.escape(correct_species_id)}\.(\d{{4}})")
    max_id = 0
    for subdir in ['Replicons', 'Genes', 'Proteins', 'RNA', 'LSTINFO']:
        folder = os.path.join(gembase_root, subdir)
        if not os.path.isdir(folder):
            continue
        for root, _, files in os.walk(folder):
            for file in files:
                match = pattern.search(file)
                if match:
                    num = int(match.group(1))
                    max_id = max(max_id, num)
    return max_id

def main():
    parser = argparse.ArgumentParser(description="Correct missassignations")
    parser.add_argument("-c", "--corrections", required=True, help="List of corrections to be performed")
    args = parser.parse_args()
    corrections_file = args.corrections
    corrections = read_corrections(corrections_file)
    for assigned_id, expected_id, base_dir in corrections:
        print(f"\n    Processing {assigned_id} -> {expected_id} ")
        if not os.path.isdir(base_dir):
            print(f"❌ Directory does not exist: {base_dir}")
            continue

        existing_ids = find_matching_genome_ids(base_dir, assigned_id)
        if not existing_ids: 
            print(f"❌ No genome files found for {assigned_id} in subfolders.")
            continue

        print(f"✅ Found genome IDs: {existing_ids}")
        last = get_last_assigned_number(base_dir, expected_id)
        print(f" Last used index for {expected_id} is {last:04d}")

        new_id_map = {}

        for old_id in existing_ids:
            last += 1
            new_id = f"{expected_id}.{last:04d}"
            new_id_map[old_id] = new_id

            print(f"{old_id} -> {new_id}")




if __name__ == "__main__":
    main()
            
            