import argparse
import re
from Bio import SeqIO

def build_key_mapping(fasta_file):
    key_to_newid = {}
    new_id_to_header = {}
    old_id_mapping = {}
    organism_mapping = {}
    replicon_type_mapping = {}

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

        # Parse ORGANISM (2nd, 3rd, 4th words of the header)
        header_parts = header.split()
        if len(header_parts) >= 4:
            organism = f"{header_parts[1]} {header_parts[2]} {header_parts[3].rstrip(',')}"
            organism_mapping[new_id] = organism
        else:
            raise ValueError(f"❌ Cannot extract organism from header: {header}")

        # Determine replicon type based on header
        if "|Complete chromosome|" in header:
            replicon_type = "complete genome"
        elif "|Complete plasmid|" in header:
            replicon_type = "complete plasmid"
        elif "|Uncultured MAG|" in header:
            replicon_type = "uncultured MAG"
        else:
            replicon_type = "draft genome"

        replicon_type_mapping[new_id] = replicon_type

    print(f"✅ Found {len(key_to_newid)} key mappings from FASTA.")
    return key_to_newid, new_id_to_header, old_id_mapping, organism_mapping, replicon_type_mapping

def replace_gbff_fields(gbff_file, key_to_newid, new_id_to_header, old_id_mapping,
                        organism_mapping, replicon_type_mapping, locus_tag_map, output_file):
    with open(gbff_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_new_id = None
        in_definition = False

        for line in infile:
            # Replace LOCUS ID with new ID
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

            # Replace DEFINITION with FASTA header
            elif line.startswith("DEFINITION") and current_new_id:
                fasta_header = new_id_to_header[current_new_id]
                outfile.write(f"DEFINITION  {fasta_header}\n")
                in_definition = True
                continue

            # Replace ACCESSION with old ID
            elif line.startswith("ACCESSION") and current_new_id:
                old_accession = old_id_mapping[current_new_id]
                outfile.write(f"ACCESSION   {old_accession}\n")
                continue

            # Replace VERSION
            elif line.startswith("VERSION") and current_new_id:
                outfile.write("VERSION     Gembases\n")
                continue

            # Replace ORGANISM with extracted one
            elif line.strip().startswith("ORGANISM") and current_new_id:
                organism = organism_mapping[current_new_id]
                outfile.write(f"  ORGANISM  {organism}\n")
                continue

            # Replace /mol_type=
            elif '/mol_type=' in line and current_new_id:
                mol_type = replicon_type_mapping[current_new_id]
                outfile.write(f'                     /mol_type="{mol_type}"\n')
                continue

            # Replace all occurrences of /locus_tag= and add OLD_locus_tag if replaced
            elif '/locus_tag=' in line and current_new_id:
                old_tags = []

                def replace_locus_tag(match):
                    old_tag = match.group(1)
                    new_tag = locus_tag_map.get(old_tag, old_tag)
                    if old_tag != new_tag:
                        old_tags.append(old_tag)
                    return f'/locus_tag="{new_tag}"'

                new_line = re.sub(r'/locus_tag\s*=\s*"([^"]+)"', replace_locus_tag, line)
                outfile.write(new_line)
                for old_tag in old_tags:
                    outfile.write(f'                     /OLD_locus_tag="{old_tag}"\n')
                continue

            # Replace /protein_id= with new locus tag
            elif '/protein_id=' in line and current_new_id:
                protein_id_match = re.search(r'/protein_id="([^"]+)"', line)
                if protein_id_match:
                    old_protein_id = protein_id_match.group(1)
                    old_locus_tag = old_protein_id.split('|')[-1]
                    new_locus_tag = locus_tag_map.get(old_locus_tag, old_locus_tag)
                    new_protein_id = f'gnl|Bakta|{new_locus_tag}'
                    line = line.replace(old_protein_id, new_protein_id)
                outfile.write(line)
                continue

            # Skip continuation lines of DEFINITION field
            elif in_definition:
                if not line.startswith(" "):
                    in_definition = False
                    outfile.write(line)
                continue

            else:
                outfile.write(line)

        print(f"✅ LOCUS, DEFINITION, ACCESSION, VERSION, ORGANISM, /mol_type, /locus_tag, and /protein_id replacements completed.")
        print(f"✅ Updated GBFF written to: {output_file}")        

def build_locus_tag_map(tsv_file):
    """
    Builds a mapping from OLDLocusTag to new LocusTag from the renamed tsv file.
    Expected columns:
    #Sequence Id	Type	Start	Stop	Strand	LocusTag	OLDLocusTag	Gene	Product	DbXrefs
    """
    locus_tag_map = {}
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            new_locus = parts[5]
            old_locus = parts[6]
            if old_locus and new_locus:
                locus_tag_map[old_locus] = new_locus
    print(f"✅ Built locus_tag mapping for {len(locus_tag_map)} entries from TSV.")
    return locus_tag_map
    
def main():
    parser = argparse.ArgumentParser(description="Replace LOCUS, DEFINITION, ACCESSION, VERSION, ORGANISM, and /mol_type in GBFF using renamed FASTA.")
    parser.add_argument("--fasta", required=True, help="Renamed FASTA file with |KEY:XYZ| and |OLD:XYZ| in headers.")
    parser.add_argument("--gbff", required=True, help="Original GBFF file.")
    parser.add_argument("--tsv", required=True, help="Renamed TSV file with LocusTag mappings.")
    parser.add_argument("--out", required=True, help="Output GBFF file with updated lines.")
    args = parser.parse_args()

    try:
        key_to_newid, new_id_to_header, old_id_mapping, organism_mapping, replicon_type_mapping = build_key_mapping(args.fasta)
        locus_tag_map = build_locus_tag_map(args.tsv)
        replace_gbff_fields(args.gbff, key_to_newid, new_id_to_header, old_id_mapping, organism_mapping, replicon_type_mapping, locus_tag_map, args.out)
    except Exception as e:
        print(f"❌ ERROR: {e}")

if __name__ == "__main__":
    main()