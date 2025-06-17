import argparse
from Bio import SeqIO
import re
import sys

def detect_chrom_or_plasmid(orig_header, new_header):
    is_complete = "[completeness=complete]" in new_header and "[topology=circular]" in new_header
    if is_complete and "chromosome" in orig_header:
        return "chromosome"
    elif is_complete and "plasmid" in orig_header:
        return "plasmid"
    else:
        return None

def parse_old_id(orig_header):
    match = re.search(r'\|(OLD:[^|]+)\|', orig_header)
    if not match:
        raise ValueError(f"Missing |OLD:...| identifier in header: {orig_header}")
    return match.group(1)

def parse_accession_species_strain(orig_header):
    """
    Extract accession, species, and strain. Strain can contain spaces.
    """
    match = re.match(r'^(?P<accession>\S+)\s+(?P<species>[A-Za-z]+ [a-z]+)\s+(?P<strain>.+?)(?:,|\||$)', orig_header)
    if match:
        return match.group("accession"), match.group("species"), match.group("strain").strip()
    else:
        raise ValueError(f"Could not parse accession/species/strain from header: {orig_header}")

def parse_mag_header(orig_header):
    """
    Specifically detect MAG format: accession MAG uncultured Genus sp. isolate Strain, ...
    """
    match = re.match(r'^(?P<accession>\S+)\s+MAG uncultured (?P<genus>[A-Za-z]+) sp\. isolate (?P<strain>\S+)', orig_header)
    if match:
        return match.group("accession"), f"{match.group('genus')} sp.", match.group("strain")
    return None

def reheader_fasta(original_fasta, bakta_fna, output_fasta):
    original_records = list(SeqIO.parse(original_fasta, "fasta"))
    bakta_records = list(SeqIO.parse(bakta_fna, "fasta"))

    if len(original_records) != len(bakta_records):
        raise ValueError("Number of sequences differs between files!")

    chrom_counter = 1
    plasmid_counter = 1
    mag_counter = 1
    draft_counter = 1

    with open(output_fasta, "w") as out_f:
        for idx, (orig_rec, bakta_rec) in enumerate(zip(original_records, bakta_records), 1):
            if str(orig_rec.seq) != str(bakta_rec.seq):
                raise ValueError(f"Sequence mismatch at sequence {idx}!")

            old_id = parse_old_id(orig_rec.description)
            key_id = bakta_rec.id  # <-- Use original bakta header here

            # Check chromosome or plasmid first
            new_header_type = detect_chrom_or_plasmid(orig_rec.description, bakta_rec.description)

            if new_header_type == "chromosome":
                accession, species, strain = parse_accession_species_strain(orig_rec.description)
                header = f"{accession}.c{chrom_counter:03d} {species} {strain} |Complete chromosome| |{old_id}| |KEY:{key_id}|"
                chrom_counter += 1

            elif new_header_type == "plasmid":
                accession, species, strain = parse_accession_species_strain(orig_rec.description)
                plasmid_name_match = re.search(r'(plasmid\s+\S+)', orig_rec.description, re.IGNORECASE)
                plasmid_name = plasmid_name_match.group(1) if plasmid_name_match else f"Plasmid_{plasmid_counter}"
                header = f"{accession}.p{plasmid_counter:03d} {species} {strain} {plasmid_name} |Complete plasmid| |{old_id}| |KEY:{key_id}|"
                plasmid_counter += 1

            else:
                mag_info = parse_mag_header(orig_rec.description)
                if mag_info:
                    accession, species, strain = mag_info
                    header = f"{accession}.m{mag_counter:03d} {species} {strain} |Uncultured MAG| |{old_id}| |KEY:{key_id}|"
                    mag_counter += 1
                else:
                    accession, species, strain = parse_accession_species_strain(orig_rec.description)
                    header = f"{accession}.d{draft_counter:03d} {species} {strain} |{old_id}| |KEY:{key_id}|"
                    draft_counter += 1

            bakta_rec.id = header
            bakta_rec.description = ""
            SeqIO.write(bakta_rec, out_f, "fasta")

    print(f"✅ New FASTA written to: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Reformat FASTA headers with structured chromosome/plasmid/MAG/draft scheme + |KEY| for bakta connection.")
    parser.add_argument("--old", required=True, help="Original FASTA file (used for header information).")
    parser.add_argument("--new", required=True, help="Bakta .fna file.")
    parser.add_argument("--out", required=True, help="Output FASTA file with new headers.")
    args = parser.parse_args()

    try:
        reheader_fasta(args.old, args.new, args.out)
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()