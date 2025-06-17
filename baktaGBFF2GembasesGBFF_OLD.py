import argparse
import re
from Bio import SeqIO

def build_key_mapping(fasta_file):
    key_to_newid = {}
    new_id_to_header = {}
    old_id_mapping = {}
    new_id_to_organism = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description

        key_match = re.search(r"\|KEY:([^\|]+)\|?", header)
        old_match = re.search(r"\|OLD:([^\|]+)\|?", header)

        if not key_match or not old_match:
            raise ValueError(f"❌ Missing |KEY:XYZ| or |OLD:XYZ| in header: {header}")

        key = key_match.group(1)
        new_id = header.split()[0]
        old_id = old_match.group(1)

        key_to_newid[key] = new_id
        new_id_to_header[new_id] = header
        old_id_mapping[new_id] = old_id

        # Extract Genus species strainID from header
        m = re.match(r"^\S+\s+([\w\-]+)\s+([\w\-]+)\s+([\w\-.]+)", header)
        if m:
            genus = m.group(1)
            species = m.group(2)
            strain = m.group(3)
            new_id_to_organism[new_id] = f"{genus} {species} {strain}"
        else:
            raise ValueError(f"❌ Could not parse Genus species strainID from header: {header}")

    print(f"✅ Found {len(key_to_newid)} key mappings from FASTA.")
    return key_to_newid, new_id_to_header, old_id_mapping, new_id_to_organism

def replace_gbff_fields(gbff_file, key_to_newid, new_id_to_header, old_id_mapping, new_id_to_organism, output_file):
    with open(gbff_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_new_id = None
        in_definition = False
        in_organism = False

        for line in infile:
            if line.startswith("LOCUS"):
                parts = line.split()
                if len(parts) >= 2:
                    old_seqid = parts[1]
                    if old_seqid in key_to_newid:
                        current_new_id = key_to_newid[old_seqid]
                        new_line = line.replace(old_seqid, current_new_id, 1)
                        outfile.write(new_line)
                        in_definition = False
                        continue
                    else:
                        raise ValueError(f"❌ LOCUS ID '{old_seqid}' not found in FASTA mapping.")

            elif line.startswith("DEFINITION") and current_new_id:
                fasta_header = new_id_to_header[current_new_id]
                outfile.write(f"DEFINITION  {fasta_header}\n")
                in_definition = True
                continue

            elif line.startswith("ACCESSION") and current_new_id:
                old_accession = old_id_mapping[current_new_id]
                outfile.write(f"ACCESSION   {old_accession}\n")
                continue

            elif line.startswith("VERSION") and current_new_id:
                outfile.write("VERSION     Gembases\n")
                continue

            elif line.startswith("  ORGANISM") and current_new_id:
                organism = new_id_to_organism[current_new_id]
                outfile.write(f"  ORGANISM  {organism}\n")
                in_organism = True
                continue

            elif in_definition:
                if not line.startswith(" "):
                    in_definition = False
                    outfile.write(line)
                # Skip continuation lines of DEFINITION
                continue

            elif in_organism:
                if not line.startswith(" "):
                    in_organism = False
                    outfile.write(line)
                # Skip taxonomy continuation lines
                continue

            else:
                outfile.write(line)

        print(f"✅ LOCUS, DEFINITION, ACCESSION, VERSION, and ORGANISM replacements completed.")
        print(f"✅ Updated GBFF written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Replace LOCUS, DEFINITION, ACCESSION, VERSION, and ORGANISM lines in GBFF using renamed FASTA.")
    parser.add_argument("--fasta", required=True, help="Renamed FASTA file with |KEY:XYZ| and |OLD:XYZ| in headers.")
    parser.add_argument("--gbff", required=True, help="Original GBFF file.")
    parser.add_argument("--out", required=True, help="Output GBFF file with updated lines.")
    args = parser.parse_args()

    try:
        key_to_newid, new_id_to_header, old_id_mapping, new_id_to_organism = build_key_mapping(args.fasta)
        replace_gbff_fields(args.gbff, key_to_newid, new_id_to_header, old_id_mapping, new_id_to_organism, args.out)
    except Exception as e:
        print(f"❌ ERROR: {e}")

if __name__ == "__main__":
    main()