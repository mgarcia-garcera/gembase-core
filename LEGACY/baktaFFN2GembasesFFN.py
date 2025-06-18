import argparse
import os
from collections import defaultdict

def parse_tsv(tsv_file, debug=False):
    features = {}
    file_prefix = None

    with open(tsv_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            seq_id, ftype, start, stop, strand, locus_tag, old_locus_tag, gene, product, dbxrefs = parts[:10]
            if not file_prefix:
                file_prefix = ".".join(locus_tag.split(".")[:3])
            features[old_locus_tag] = {
                'SequenceId': seq_id,
                'Type': ftype,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'LocusTag': locus_tag,
                'OLDLocusTag': old_locus_tag,
                'Gene': gene,
                'Product': product,
                'DbXrefs': dbxrefs,
            }
            if debug:
                print(f"[DEBUG] Parsed feature: {features[old_locus_tag]}")
    return features, file_prefix

def generate_new_gene_id(locus_tag, count):
    parts = locus_tag.split(".")
    if len(parts) < 3:
        return f"unknown_{count:05d}"
    aabb = parts[0][:4].lower()
    XXX = parts[1]
    YYYY = parts[2]
    ZZZZZ = f"{count:05d}"
    return f"{aabb}{XXX}.{YYYY}_{ZZZZZ}"

def build_fasta_header(feat, count_by_type):
    gene = feat['Gene']
    if not gene:
        count_by_type[feat['Type']] += 1
        gene = generate_new_gene_id(feat['LocusTag'], count_by_type[feat['Type']])
    return " ".join([
        feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
        gene, feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
    ])

def parse_fasta_headers(fasta_file):
    headers = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]  # first word
                headers.add(header)
    return headers

def check_headers(features, fasta_headers, fasta_name, debug=False):
    old_locus_tags = set(features.keys())
    found = old_locus_tags & fasta_headers
    missing = old_locus_tags - fasta_headers

    print(f"\n[INFO] Checking headers for {fasta_name}:")
    print(f"  - Total OLDLocusTags in TSV: {len(old_locus_tags)}")
    print(f"  - Found in {fasta_name}      : {len(found)}")
    print(f"  - Missing in {fasta_name}    : {len(missing)}")

    if missing:
        print(f"  [WARNING] Missing OLDLocusTags in {fasta_name}:")
        for tag in sorted(missing):
            print(f"    {tag}")

    if debug and found:
        print(f"[DEBUG] Example of matched locus tags: {list(found)[:5]}")

def write_fasta(outfile, header, seq):
    with open(outfile, 'a') as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def process_fasta(input_file, features, file_prefix, suffix_map, out_dirs, count_by_type, debug=False):
    with open(input_file) as handle:
        seq_lines = []
        header = None

        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    old_locus = header[1:].split()[0]
                    if old_locus in features:
                        feat = features[old_locus]
                        fasta_header = build_fasta_header(feat, count_by_type)
                        suffix = suffix_map.get(feat['Type'].lower(), '.misc')
                        out_dir = out_dirs['nuc'] if suffix == '.nuc' else out_dirs['rna']
                        out_file = os.path.join(out_dir, f"{file_prefix}{suffix}")
                        write_fasta(out_file, fasta_header, "".join(seq_lines))
                        #if debug:
                        #    print(f"[DEBUG] Wrote {old_locus} to {out_file}")
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            old_locus = header[1:].split()[0]
            if old_locus in features:
                feat = features[old_locus]
                fasta_header = build_fasta_header(feat, count_by_type)
                suffix = suffix_map.get(feat['Type'].lower(), '.misc')
                out_dir = out_dirs['nuc'] if suffix == '.nuc' else out_dirs['rna']
                out_file = os.path.join(out_dir, f"{file_prefix}{suffix}")
                write_fasta(out_file, fasta_header, "".join(seq_lines))
                #if debug:
                #    print(f"[DEBUG] Wrote {old_locus} to {out_file}")

def process_faa(input_file, features, file_prefix, protein_dir, count_by_type, debug=False):
    outfile = os.path.join(protein_dir, f"{file_prefix}.faa")
    open(outfile, 'w').close()
    with open(input_file) as handle:
        seq_lines = []
        header = None

        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    old_locus = header[1:].split()[0]
                    if old_locus in features:
                        feat = features[old_locus]
                        fasta_header = build_fasta_header(feat, count_by_type)
                        write_fasta(outfile, fasta_header, "".join(seq_lines))
                        #if debug:
                        #    print(f"[DEBUG] Wrote {old_locus} to {outfile}")
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            old_locus = header[1:].split()[0]
            if old_locus in features:
                feat = features[old_locus]
                fasta_header = build_fasta_header(feat, count_by_type)
                write_fasta(outfile, fasta_header, "".join(seq_lines))
                #if debug:
                #    print(f"[DEBUG] Wrote {old_locus} to {outfile}")

def main(tsv_file, ffn_file, faa_file, out_dir, debug=False):
    features, file_prefix = parse_tsv(tsv_file, debug=debug)

    print("\n========== HEADER CONSISTENCY CHECK ==========")
    ffn_headers = parse_fasta_headers(ffn_file)
    check_headers(features, ffn_headers, os.path.basename(ffn_file), debug=debug)
    faa_headers = parse_fasta_headers(faa_file)
    check_headers(features, faa_headers, os.path.basename(faa_file), debug=debug)
    print("==============================================\n")

    nuc_dir = os.path.join(out_dir, "nucleotide")
    rna_dir = os.path.join(out_dir, "RNA")
    protein_dir = os.path.join(out_dir, "protein")
    os.makedirs(nuc_dir, exist_ok=True)
    os.makedirs(rna_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)

    suffix_map = {
        'cds': '.nuc',
        'rrna': '.rrna',
        'trna': '.trna',
        'tmrna': '.tmrna',
        'ncrna': '.ncrna',
        'ncrna-region': '.ncrna'
    }

    out_dirs = {'nuc': nuc_dir, 'rna': rna_dir}
    count_by_type = defaultdict(int)

    print(f"[INFO] Processing FFN file...")
    process_fasta(ffn_file, features, file_prefix, suffix_map, out_dirs, count_by_type, debug=debug)

    print(f"[INFO] Processing FAA file...")
    process_faa(faa_file, features, file_prefix, protein_dir, count_by_type, debug=debug)

    print(f"[INFO] Completed. Files written to:")
    print(f"       {nuc_dir}  (nucleotide)")
    print(f"       {rna_dir}  (RNA)")
    print(f"       {protein_dir}  (protein)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rewrite FASTA headers from TSV, with header consistency checking.")
    parser.add_argument("--tsv", required=True, help="TSV file with features")
    parser.add_argument("--ffn", required=True, help="FFN nucleotide FASTA file")
    parser.add_argument("--faa", required=True, help="FAA protein FASTA file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()

    main(args.tsv, args.ffn, args.faa, args.outdir, debug=args.debug)

