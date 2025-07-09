import os
import re
import argparse
import shutil
import gzip

def parse_tsv(tsv_path):
    species_map = {}
    with open(tsv_path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                wrong_id, species_name, correct_id, location = line.strip().split("\t")
                key = (species_name.strip(), correct_id.strip(), location.strip())
                if key not in species_map:
                    species_map[key] = []
                species_map[key].append(wrong_id.strip())
    return species_map

def get_last_assigned_index(base_path, correct_id):
    max_index = 0
    pattern = re.compile(rf"{re.escape(correct_id)}\.(\d{{4}})")
    for root, _, files in os.walk(base_path):
        for fname in files:
            match = pattern.search(fname)
            if match:
                idx = int(match.group(1))
                max_index = max(max_index, idx)
            elif fname.endswith(".gz"):
                try:
                    with gzip.open(os.path.join(root, fname), 'rt', encoding='utf-8') as f:
                        for line in f:
                            for m in pattern.finditer(line):
                                idx = int(m.group(1))
                                max_index = max(max_index, idx)
                except Exception:
                    continue
    return max_index

def update_file_content(path, old_id, new_id, is_gz=False):
    try:
        if is_gz:
            with gzip.open(path, 'rt', encoding='utf-8') as f:
                content = f.read()
            content = re.sub(rf"{re.escape(old_id)}", new_id, content)
            with gzip.open(path, 'wt', encoding='utf-8') as f:
                f.write(content)
        else:
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()
            content = re.sub(rf"{re.escape(old_id)}", new_id, content)
            with open(path, 'w', encoding='utf-8') as f:
                f.write(content)
    except Exception as e:
        print(f"Failed to update file {path}: {e}")

def process_species_group(species_name, correct_id, location, wrong_ids):
    print(f"\nProcessing species: {species_name}")
    print(f"  Correct ID: {correct_id}")
    print(f"  Directory: {location}")
    print(f"  Wrong IDs: {wrong_ids}")
    
    counter = get_last_assigned_index(location, correct_id)
    print(f"  Starting from index: {counter:04d}")

    pattern = re.compile(r"(" + "|".join(re.escape(wid) for wid in wrong_ids) + r")\.(\d{4})")
    
    for root, _, files in os.walk(location):
        for fname in files:
            match = pattern.search(fname)
            if not match:
                continue
            old_prefix = match.group(1)
            old_index = match.group(2)
            old_id = f"{old_prefix}.{old_index}"
            counter += 1
            new_id = f"{correct_id}.{counter:04d}"

            old_path = os.path.join(root, fname)
            new_fname = pattern.sub(new_id, fname)
            new_path = os.path.join(root, new_fname)

            print(f"  Renaming {fname} → {new_fname}")
            shutil.move(old_path, new_path)

            is_gz = new_path.endswith(".gz")
            update_file_content(new_path, old_id, new_id, is_gz=is_gz)

def main():
    parser = argparse.ArgumentParser(description="Fix genome IDs and reassign them properly.")
    parser.add_argument("-c","--correction", help="Path to the correction TSV file")
    args = parser.parse_args()

    species_map = parse_tsv(args.correction)

    for (species_name, correct_id, location), wrong_ids in species_map.items():
        process_species_group(species_name, correct_id, location, wrong_ids)

if __name__ == "__main__":
    main()