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

def parse_mag_header(orig_header):
    match = re.match(r'^(?P<accession>\S+)\s+MAG uncultured (?P<species>.+?) sp\. isolate .*?_bin\.(?P<bin>\d+).*?\|OLD:[^|]+\|', orig_header)
    if match:
        d = match.groupdict()
        d['bin'] = int(d['bin'])
        return d

    match = re.match(r'^(?P<accession>\S+)\s+MAG uncultured (?P<species>.+?) sp\. isolate (?P<isolate>\S+).*?\|OLD:[^|]+\|', orig_header)
    if match:
        d = match.groupdict()
        d['bin'] = None
        return d

    return None

def parse_default_header(orig_header):
    match = re.match(r'^(?P<accession>\S+)\s+(?P<species>.+?) strain (?P<strain>\S+).*NODE_\d+_length_(?P<length>\d+)_cov_(?P<cov>[\d\.]+).*?\|OLD:[^|]+\|', orig_header)
    if match:
        return match.groupdict()

    match = re.match(r'^(?P<accession>\S+)\s+(?P<species>.+?) (?P<strain>\S+)\s+(supercontig|supercont\d*|scaffold\d*|contig\d*|supercont\d+\.\d+)[^\|]*?\|OLD:[^|]+\|', orig_header)
    if match:
        gd = match.groupdict()
        gd['length'] = None
        gd['cov'] = None
        return gd

    match = re.match(r'^(?P<accession>\S+)\s+(?P<species>.+?) (?P<strain>\S+).*whole genome shotgun sequence.*?\|OLD:[^|]+\|', orig_header)
    if match:
        gd = match.groupdict()
        gd['length'] = None
        gd['cov'] = None
        return gd

    return None

def extract_species_and_strain(orig_header):
    match = re.match(r'^\S+\s+(?P<species>.+?) strain (?P<strain>\S+)', orig_header)
    if match:
        return match.group('species'), match.group('strain')
    else:
        return "Unknown_species", "Unknown_strain"

def reheader_fasta(original_fasta, bakta_fna, output_fasta):
    original_records = list(SeqIO.parse(original_fasta, "fasta"))
    bakta_records = list(SeqIO.parse(bakta_fna, "fasta"))

    if len(original_records) != len(bakta_records):
        raise ValueError("Number of sequences differs between files!")

    chrom_counter = 1
    plasmid_counter = 1
    bin_counter = 1

    with open(output_fasta, "w") as out_f:
        for idx, (orig_rec, bakta_rec) in enumerate(zip(original_records, bakta_records), 1):
            if str(orig_rec.seq) != str(bakta_rec.seq):
                raise ValueError(f"Sequence mismatch at sequence {idx}!")

            old_id = parse_old_id(orig_rec.description)
            new_header_type = detect_chrom_or_plasmid(orig_rec.description, bakta_rec.description)
            accession = orig_rec.id

            if new_header_type == "chromosome":
                species, strain = extract_species_and_strain(orig_rec.description)
                new_header = f"CHROMOSOME_{chrom_counter} ({accession}) {species} {strain} |Complete chromosome| |{old_id}|"
                chrom_counter += 1

            elif new_header_type == "plasmid":
                species, strain = extract_species_and_strain(orig_rec.description)
                new_header = f"PLASMID_{plasmid_counter} ({accession}) {species} {strain} |Complete plasmid| |{old_id}|"
                plasmid_counter += 1

            else:
                mag_info = parse_mag_header(orig_rec.description)
                if mag_info:
                    if mag_info['bin'] is None:
                        assigned_bin = bin_counter
                        bin_counter += 1
                    else:
                        assigned_bin = mag_info['bin']
                    new_header = f"Bin{assigned_bin} ({mag_info['accession']}) {mag_info['species']} sp. |Uncultured MAG| |{old_id}|"
                else:
                    default_info = parse_default_header(orig_rec.description)
                    if default_info:
                        if default_info.get('length'):
                            new_header = f"Contig{idx} ({default_info['accession']}) {default_info['species']} {default_info['strain']} |length:{default_info['length']}| |cov:{default_info['cov']}| |{old_id}|"
                        else:
                            new_header = f"Contig{idx} ({default_info['accession']}) {default_info['species']} {default_info['strain']} |{old_id}|"
                    else:
                        raise ValueError(f"Could not parse header for contig {idx}: {orig_rec.description}")

            bakta_rec.id = new_header
            bakta_rec.description = ""
            SeqIO.write(bakta_rec, out_f, "fasta")

    print(f"✅ New FASTA written to: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Standardize FASTA headers for chromosomes, plasmids, MAGs, or drafts. Enforces |OLD:...| presence.")
    parser.add_argument("--old", required=True, help="Original FASTA file.")
    parser.add_argument("--new", required=True, help="Bakta .fna file.")
    parser.add_argument("--out", required=True, help="Output FASTA file with standardized headers.")
    args = parser.parse_args()

    try:
        reheader_fasta(args.old, args.new, args.out)
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()