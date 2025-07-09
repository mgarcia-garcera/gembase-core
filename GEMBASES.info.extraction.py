from pathlib import Path
import csv
import argparse


parser = argparse.ArgumentParser(description="Cluster species FASTAs with MMseqs2")
parser.add_argument("-i", "--input_folder", required=True, help="input folder")
parser.add_argument("-o", "--output", required=True, help="output file")
args = parser.parse_args()

lstinfo_dir = Path(args.input_folder)
output_csv = args.output

rows = []

for inf_file in lstinfo_dir.glob("*.inf"):
    data = {}
    with inf_file.open() as f:
        for line in f:
            if ":" not in line:
                continue
            key, value = line.strip().split(":", 1)
            data[key.strip()] = value.strip()

    # Skip if required keys are missing
    if not all(k in data for k in ("GEBMASES_Id", "NAME", "species")):
        print(f"⚠️ Skipping {inf_file.name}: missing one of required keys")
        continue

    full_id = data["GEBMASES_Id"]
    strain_name = data["NAME"]
    species_name = data["species"]
    species_id = ".".join(full_id.split(".")[:2])

    rows.append([species_id, species_name, full_id, strain_name])

# Write to tab-separated file
with open(output_csv, "a", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(rows)

print(f"✅ Summary written to {output_csv}")