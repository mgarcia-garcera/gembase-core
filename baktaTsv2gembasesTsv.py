import argparse
import re
import csv
from collections import defaultdict
from Bio import SeqIO

def build_key_mapping(fasta_file):
    mapping = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        key_match = re.search(r"\|KEY:([^\|]+)\|?", header)
        if key_match:
            key = key_match.group(1)
            new_id = header.split()[0]
            mapping[key] = new_id
        else:
            raise ValueError(f"❌ Missing |KEY:XYZ| in header: {header}")
    print(f"✅ Found {len(mapping)} key mappings from FASTA file.")
    return mapping

def reformat_tsv(tsv_file, mapping, output_file):
    with open(tsv_file, 'r', newline='') as in_f:
        lines = in_f.readlines()

    header_line_index = None
    for i, line in enumerate(lines):
        if line.startswith("#Sequence Id"):
            header_line_index = i
            break

    if header_line_index is None:
        raise ValueError("❌ Could not find header line starting with '#Sequence Id' in TSV.")

    header_line = lines[header_line_index]
    data_lines = lines[header_line_index + 1:]

    reader = csv.DictReader([header_line] + data_lines, delimiter='\t')

    # New fieldnames: insert 'LocusTag' after 'Strand' and rename 'Locus Tag' → 'OLDLocusTag'
    fieldnames = []
    for fn in reader.fieldnames:
        fieldnames.append(fn)
        if fn == "Strand":
            fieldnames.append("LocusTag")
    fieldnames = [f if f != "Locus Tag" else "OLDLocusTag" for f in fieldnames]

    locus_counters = defaultdict(int)

    with open(output_file, 'w', newline='') as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        updated_rows = 0

        for row in reader:
            old_seqid = row["#Sequence Id"]
            if old_seqid not in mapping:
                raise ValueError(f"❌ SeqID '{old_seqid}' from TSV not found in FASTA mapping!")

            new_seqid = mapping[old_seqid]
            row["#Sequence Id"] = new_seqid

            # Replace strand symbols
            strand = row["Strand"]
            if strand == "+":
                row["Strand"] = "D"
            elif strand == "-":
                row["Strand"] = "C"
            else:
                raise ValueError(f"❌ Unknown strand value '{strand}' in row.")

            # Generate new LocusTag
            locus_counters[new_seqid] += 1
            locus_number = f"{locus_counters[new_seqid]:05}0"
            row["LocusTag"] = f"{new_seqid}_{locus_number}"

            # Rename original 'Locus Tag' → 'OLDLocusTag'
            row["OLDLocusTag"] = row.get("Locus Tag", "")

            # Clean up to prevent dict mismatch errors
            if "Locus Tag" in row:
                del row["Locus Tag"]

            writer.writerow(row)
            updated_rows += 1

    print(f"✅ Reformatted TSV written to: {output_file} ({updated_rows} rows updated)")

def main():
    parser = argparse.ArgumentParser(description="Reformat Bakta TSV: update SeqIDs, replace strand +/-, add LocusTag, rename old Locus Tag to OLDLocusTag.")
    parser.add_argument("--fasta", required=True, help="Renamed FASTA file with |KEY:XYZ| in headers.")
    parser.add_argument("--tsv", required=True, help="Bakta output TSV file.")
    parser.add_argument("--out", required=True, help="Output TSV with updated SeqIDs and LocusTags.")
    args = parser.parse_args()

    try:
        mapping = build_key_mapping(args.fasta)
        reformat_tsv(args.tsv, mapping, args.out)
    except Exception as e:
        print(f"❌ ERROR: {e}")

if __name__ == "__main__":
    main()