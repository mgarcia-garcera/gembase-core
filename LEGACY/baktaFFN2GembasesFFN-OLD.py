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

def write_fasta(outfile, header, seq):
    with open(outfile, 'a') as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def get_outfile_path(folder, file_prefix, suffix):
    return os.path.join(folder, f"{file_prefix}{suffix}")

def process_fasta(input_file, features, file_prefix, suffix_map, nuc_dir, rna_dir, out_files, no_gene_count, debug=False):
    with open(input_file) as handle:
        seq_lines = []
        header = None

        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    old_locus = header.split()[0].lstrip(">")
                    if old_locus not in features:
                        if debug:
                            print(f"[DEBUG] OLDLocusTag {old_locus} not found in TSV. Skipping.")
                    else:
                        feat = features[old_locus]
                        gene = feat['Gene']
                        if not gene:
                            no_gene_count[feat['Type']] += 1
                            gene = generate_new_gene_id(feat['LocusTag'], no_gene_count[feat['Type']])

                        new_header = " ".join([
                            feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
                            gene, feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
                        ])
                        suffix = suffix_map.get(feat['Type'].lower(), '.misc')
                        if suffix == '.nuc':
                            out_file = get_outfile_path(nuc_dir, file_prefix, suffix)
                        else:
                            out_file = get_outfile_path(rna_dir, file_prefix, suffix)

                        if out_file not in out_files:
                            open(out_file, 'w').close()
                            out_files[out_file] = True

                        write_fasta(out_file, new_header, "".join(seq_lines))
                        if debug:
                            print(f"[DEBUG] Wrote {old_locus} to {out_file}")

                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            old_locus = header.split()[0].lstrip(">")
            if old_locus in features:
                feat = features[old_locus]
                gene = feat['Gene']
                if not gene:
                    no_gene_count[feat['Type']] += 1
                    gene = generate_new_gene_id(feat['LocusTag'], no_gene_count[feat['Type']])
                new_header = " ".join([
                    feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
                    gene, feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
                ])
                suffix = suffix_map.get(feat['Type'].lower(), '.misc')
                if suffix == '.nuc':
                    out_file = get_outfile_path(nuc_dir, file_prefix, suffix)
                else:
                    out_file = get_outfile_path(rna_dir, file_prefix, suffix)

                if out_file not in out_files:
                    open(out_file, 'w').close()
                    out_files[out_file] = True

                write_fasta(out_file, new_header, "".join(seq_lines))
                if debug:
                    print(f"[DEBUG] Wrote {old_locus} to {out_file}")

def process_faa(input_file, features, file_prefix, protein_dir, debug=False):
    outfile = os.path.join(protein_dir, f"{file_prefix}.faa")
    open(outfile, 'w').close()
    with open(input_file) as handle:
        seq_lines = []
        header = None

        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    old_locus = header.split()[0].lstrip(">")
                    if old_locus not in features:
                        if debug:
                            print(f"[DEBUG] OLDLocusTag {old_locus} not found in TSV. Skipping.")
                    else:
                        feat = features[old_locus]
                        gene = feat['Gene']
                        new_header = " ".join([
                            feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
                            gene if gene else "NA", feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
                        ])
                        write_fasta(outfile, new_header, "".join(seq_lines))
                        if debug:
                            print(f"[DEBUG] Wrote {old_locus} to {outfile}")

                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            old_locus = header.split()[0].lstrip(">")
            if old_locus in features:
                feat = features[old_locus]
                gene = feat['Gene']
                new_header = " ".join([
                    feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
                    gene if gene else "NA", feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
                ])
                write_fasta(outfile, new_header, "".join(seq_lines))
                if debug:
                    print(f"[DEBUG] Wrote {old_locus} to {outfile}")

def main(tsv_file, ffn_file, faa_file, out_dir, debug=False):
    features, file_prefix = parse_tsv(tsv_file, debug=debug)

    suffix_map = {
        'cds': '.nuc',
        'rrna': '.rrna',
        'trna': '.ttrna',
        'tmrna': '.tmrna',
        'ncrna': '.ncrna',
        'ncrna-region': '.ncrna'
    }

    nuc_dir = os.path.join(out_dir, "Genes")
    rna_dir = os.path.join(out_dir, "RNA")
    protein_dir = os.path.join(out_dir, "Proteins")
    os.makedirs(nuc_dir, exist_ok=True)
    os.makedirs(rna_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)

    out_files = {}
    no_gene_count = defaultdict(int)

    print(f"[INFO] Processing FFN file...")
    process_fasta(ffn_file, features, file_prefix, suffix_map, nuc_dir, rna_dir, out_files, no_gene_count, debug=debug)
    
    print(f"[INFO] Processing FAA file...")
    process_faa(faa_file, features, file_prefix, protein_dir, debug=debug)

    print(f"[INFO] Completed. Files written to:")
    print(f"       {nuc_dir}  (for .nuc)")
    print(f"       {rna_dir}  (for .rrna, .ttrna, etc.)")
    print(f"       {protein_dir}  (for protein FASTA)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rewrite FFN & FAA FASTA headers from TSV metadata and organize by type.")
    parser.add_argument("--tsv", required=True, help="Renamed TSV input file")
    parser.add_argument("--ffn", required=True, help="FFN nucleotide FASTA file")
    parser.add_argument("--faa", required=True, help="FAA protein FASTA file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory where subfolders will be created")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")

    args = parser.parse_args()
    
    #main(args.tsv, args.ffn, args.faa, args.outdir, debug=args.debug)rs will be created")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()

    main(args.tsv, args.ffn, args.faa, args.outdir, debug=args.debug)