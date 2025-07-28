import struct
import Bio
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import os
import re
import subprocess
import shutil
import sys
from collections import defaultdict
from collections import Counter
import pathlib
from pathlib import Path
import gzip

def extract_definition_keywords(gbff_path):
    """
    Extracts unique keywords from DEFINITION lines for each LOCUS in a GBFF.

    Returns:
        dict: {LOCUS_ID (as in |OLD:XYZ|) → unique keyword string}
    """
    open_func = gzip.open if gbff_path.endswith('.gz') else open

    locus_definitions = {}  # LOCUS → DEFINITION
    token_lists = {}        # LOCUS → set(tokens)

    with open_func(gbff_path, 'rt') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            locus_id = record.id.strip()
            definition = record.description.strip()
            tokens = set(re.findall(r'\b\w+\b', definition.lower()))  # simple word tokenizer

            locus_definitions[locus_id] = definition
            token_lists[locus_id] = tokens

    # Compute common keywords across all definitions
    all_token_sets = list(token_lists.values())
    if not all_token_sets:
        return {}

    common_keywords = set.intersection(*all_token_sets)

    # Compute unique term(s) per LOCUS by subtracting common ones
    locus_unique = {}
    for locus, tokens in token_lists.items():
        unique_tokens = sorted(tokens - common_keywords)
        locus_unique[locus] = " ".join(unique_tokens) if unique_tokens else ""

    return locus_unique

def extract_fna_from_gbff(gbff_path):
    """
    Extracts genomic sequences from a GBFF or GBFF.GZ file and writes them to a temporary FASTA (fna) file.
    
    Parameters:
        gbff_path (str): Path to the GBFF file (can be .gbff or .gbff.gz).
    
    Returns:
        str: Path to the temporary fna (FASTA) file.
    """
    try:
        # Open depending on whether file is gzipped
        if gbff_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'  # read text mode
        else:
            open_func = open
            mode = 'r'

        temp_fna = tempfile.NamedTemporaryFile(delete=False, suffix=".fna", mode='w')

        with open_func(gbff_path, mode) as gbff_handle:
            records = SeqIO.parse(gbff_handle, "genbank")
            count = SeqIO.write(records, temp_fna, "fasta")
        
        temp_fna.close()
        print(f"Extracted {count} sequences to: {temp_fna.name}")
        return temp_fna.name

    except Exception as e:
        print(f"Error extracting FNA from GBFF: {e}")
        return None

def detect_file_format(file_path):
    """
    Detects whether a file is in FASTA or GBFF (GenBank Flat File Format).
    
    Parameters:
        file_path (str): Path to the file.
    
    Returns:
        str: 'FASTA', 'GBFF', or 'Unknown'
    """
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    return 'FASTA'
                elif line.startswith('LOCUS'):
                    return 'GBFF'
                else:
                    break  # Only look at the beginning
    except Exception as e:
        print(f"Error reading file: {e}")
        return 'Unknown'
    
    # Fallback to file extension if content didn't help
    if file_path.endswith(('.fasta', '.fa')):
        return 'FASTA'
    elif file_path.endswith(('.gbff', '.gb')):
        return 'GBFF'

    return 'Unknown'

def ensure_directories(base_dir, subdirs):
    for subdir in subdirs:
        full_path = os.path.join(base_dir, subdir)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
            print(f"Created directory: {full_path}")
        else:
            print(f"Directory already exists: {full_path}")

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

def parse_accession_species_strain_draft(orig_header):
    """
    Extract accession, species, and strain. Strain can contain spaces.
    """
    match = re.match(r'^(?P<accession>\S+)\s+(?P<species>[A-Za-z]+ [a-z]+)\s+(?P<strain>.+?)(?:,|\||$)', orig_header)
    if match:
        return match.group("accession"), match.group("species"), match.group("strain").strip()
    else:
        raise ValueError(f"Could not parse accession/species/strain from header in this draft: {orig_header}")

def parse_mag_header(orig_header):
    """
    Detects if 'MAG' exists in the header.
    If found, extract accession, Genus species, and strain.

    Returns:
        accession, species, strain  |  or  None if no MAG detected.
    """
    match = re.match(r'^(?P<accession>\S+)\s+(?P<desc>.+)', orig_header)
    if not match:
        return None

    accession = match.group("accession")
    desc = match.group("desc")

    if "MAG" not in desc:
        return None

    # Extract first occurrence of Genus species (e.g., "Abiotrophia defectiva")
    species_match = re.search(r'([A-Z][a-z]+)\s+(sp\.|[a-z]+)', desc)
    if not species_match:
        return None

    genus = species_match.group(1)
    species = species_match.group(2)
    species_full = f"{genus} {species}"

    # Everything AFTER that match is treated as strain (optional cleanup possible)
    strain = desc[species_match.end():].strip() or "no_strain"

    return accession, species_full, strain

def reheader_fastaGBFF(input_fasta, output_fasta):
    """
    Reheaders a single FASTA file based on replicon type and structured identifier format.
    Input format:
    >ABUT.001.0002 Genus species strain ContigID,... |OLD:XYZ|
    
    Output format:
    >ABUT.001.0002.d01 Genus species strain ContigID |OLD:XYZ|
    """

    draft_counter = 1
    chrom_counter = 1
    plasmid_counter = 1
    mag_counter = 1

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            orig_header = record.description

            # Extract old accession
            old_id = parse_old_id(orig_header)

            # Parse species/strain/etc
            accession, species, strain = parse_accession_species_strain_draft(orig_header)

            # Remove |OLD:...| and leading accession from description
            clean_desc = re.sub(r'\|OLD:[^|]+\|?', '', orig_header).strip()
            clean_desc = clean_desc.replace(',', '')

            # Remove the accession (first token)
            clean_desc = clean_desc.split(maxsplit=1)
            clean_desc = clean_desc[1] if len(clean_desc) > 1 else ""

            # Determine genome type suffix
            genome_type = detect_chrom_or_plasmid(orig_header, orig_header)
            if genome_type == "chromosome":
                suffix = f".c{chrom_counter:02d}"
                chrom_counter += 1
            elif genome_type == "plasmid":
                suffix = f".p{plasmid_counter:02d}"
                plasmid_counter += 1
            elif parse_mag_header(orig_header):
                suffix = f".m{mag_counter:02d}"
                mag_counter += 1
            else:
                suffix = f".d{draft_counter:02d}"
                draft_counter += 1

            # Construct new header
            new_id = f"{accession}{suffix}"
            new_header = f">{new_id} {clean_desc.strip()} |OLD:{old_id}|"

            outfile.write(new_header + '\n')
            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                outfile.write(seq[i:i+60] + '\n')

    print(f"✅ Reheadered FASTA written to: {output_fasta}")

def extract_features_from_gbff(gbff_path, identifier, output_tsv, replicon_mapping):
    """
    Extracts annotated features from a GBFF or GBFF.GZ file and writes them to a TSV.

    Columns:
    - LocusTag (standardized, full format)
    - Type (e.g., CDS, rRNA)
    - Start
    - Stop
    - Strand (D or C)
    - Replicon (e.g., ABUT.001.0002.d01)
    - OLDLocusTag (from original annotation)
    - Product (description)
    """

    is_gzipped = gbff_path.endswith(".gz")
    open_func = gzip.open if is_gzipped else open

    feature_counters = {}  # ✅ Per-replicon counters
    rows = []

    with open_func(gbff_path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            old_locus_id = record.id.split('.')[0]  # e.g., NZ_NIGF01000001
            if old_locus_id not in replicon_mapping:
                raise ValueError(f"❌ Replicon ID '{old_locus_id}' not found in mapping!")

            replicon = replicon_mapping[old_locus_id]  # ✅ e.g., ABUT.001.0002.d01

            # ✅ Extract suffix (e.g., .d01)
            suffix_match = re.search(r'\.([cpdm]\d{2,3})$', replicon)
            if not suffix_match:
                raise ValueError(f"❌ Cannot extract [cpdm]ZZ from replicon name: {replicon}")
            suffix = suffix_match.group(1)

            # ✅ Initialize per-replicon feature count
            if suffix not in feature_counters:
                feature_counters[suffix] = 1

            for feature in record.features:
                if feature.type.lower() not in ['cds', 'rrna', 'trna']:
                    continue

                # Coordinates
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = feature.location.strand
                strand_symbol = 'D' if strand == 1 else 'C'

                # Original locus tag
                old_tag = feature.qualifiers.get('locus_tag', ["-"])[0]

                # Product
                product = feature.qualifiers.get('product', ["-"])[0]

                # ✅ Generate locus tag using new format
                feature_number = f"{feature_counters[suffix]:06}"
                new_tag = f"{identifier}.{suffix}_{feature_number}"

                # ✅ Increment counter
                feature_counters[suffix] += 1

                # Type
                ftype = feature.type

                rows.append([
                    new_tag, ftype, start, end, strand_symbol, replicon, old_tag, product
                ])

    # Write TSV
    with open(output_tsv, "w", newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow([
            "#LocusTag", "Type", "Start", "Stop", "Strand",
            "Replicon", "OLDLocusTag", "Product"
        ])
        writer.writerows(rows)

    print(f"✅ Extracted {sum(feature_counters.values())} features to: {output_tsv}")
   
def build_old_to_new_replicon_mapping(fasta_file):
    """
    Builds a dictionary mapping from OLD accession ID to new replicon ID from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file with headers like:
                          >ABUT.001.0002.d01 ... |OLD:NZ_ABC123456|

    Returns:
        dict: { OLD_ID → NEW_ID }
    """
    mapping = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description

        # Extract OLD ID
        old_match = re.search(r'\|OLD:([^|]+)\|?', header)
        if not old_match:
            raise ValueError(f"Missing |OLD:...| in header: {header}")
        old_id = old_match.group(1)

        # Extract new ID (first word in header)
        new_id = header.split()[0]

        mapping[old_id] = new_id

    print(f"✅ Built mapping for {len(mapping)} replicons.")
    return mapping

def build_old_to_new_replicon_mapping_new(fasta_file):
    """
    Builds a dictionary mapping from OLD accession ID to new replicon ID from a FASTA file.

    Returns:
        dict: { OLD_ID → NEW_ID }
    """
    mapping = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        match = re.search(r'\|OLD:([^\|]+)\|?', header)
        if not match:
            raise ValueError(f"❌ Missing |OLD:...| in header: {header}")
        old_id = match.group(1).strip()  # ← Keep the .1 suffix if present
        new_id = header.split()[0]
        mapping[old_id] = new_id

    print(f"✅ Built mapping for {len(mapping)} replicons.")
    return mapping

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
                    accession, species, strain = parse_accession_species_strain_draft(orig_rec.description)
                    header = f"{accession}.d{draft_counter:03d} {species} {strain} |{old_id}| |KEY:{key_id}|"
                    draft_counter += 1

            bakta_rec.id = header
            bakta_rec.description = ""
            SeqIO.write(bakta_rec, out_f, "fasta")

    print(f"✅ New FASTA written to: {output_fasta}")

def br(gi, file_path):
    """
    Reads a 4-byte big-endian integer from a binary file at a position determined by `gi`.

    Parameters:
        gi (int): Index used to calculate byte position (gi * 4).
        file_path (str): Path to the binary file.

    Returns:
        int: The 4-byte integer value read from the file.
    
    Raises:
        ValueError: If fewer than 4 bytes are read (e.g. end of file).
        IOError: If the file cannot be opened or read.
    """
    pos = gi * 4
    with open(file_path, 'rb') as f:
        f.seek(pos)
        buffer = f.read(4)
        if len(buffer) < 4:
            raise ValueError(f"Unexpected end of file or insufficient data at position {pos}.")
        return struct.unpack('>I', buffer)[0]

def extract_genus(species):
    """
    Extracts the genus from a species name.
    """
    return species.split()[0] if species else None

def taxline(species, genera, taxnames, nodes, phylum_info, class_info, order_info, family_info):
    """
    Constructs a taxonomy line from a species name.

    Parameters:
        species (str): Full species name (e.g., "Homo sapiens").
        genera (dict): Map from genus to taxid.
        taxnames (dict): Map from taxid to name.
        nodes (str): Path to the binary file used in br().
        phylum_info, class_info, order_info, family_info (str): Delimited strings (e.g. "|123|456|") containing taxids.

    Returns:
        str: A taxonomy line like "Phylum, Class, Order, Family, Genus, Species".
    """
    genus = extract_genus(species)
    taxid = genera.get(genus, 0)
    taxinf = {}

    while taxid and taxid not in (0, 1):
        taxid = br(taxid, nodes)
        tag = f"|{taxid}|"
        if tag in phylum_info:
            taxinf['phylum'] = taxnames.get(taxid, "")
        elif tag in class_info:
            taxinf['class'] = taxnames.get(taxid, "")
        elif tag in order_info:
            taxinf['order'] = taxnames.get(taxid, "")
        elif tag in family_info:
            taxinf['family'] = taxnames.get(taxid, "")

    return ', '.join([
        taxinf.get('phylum', ''),
        taxinf.get('class', ''),
        taxinf.get('order', ''),
        taxinf.get('family', ''),
        genus or '',
        species or ''
    ])


def seqtrans(dna_seq,code):
    return dna_seq.translate(table=code)

def detect_file_format(file_path):
    """
    Detects whether a file is in FASTA or GBFF (GenBank Flat File Format).
    
    Parameters:
        file_path (str): Path to the file.
    
    Returns:
        str: 'FASTA', 'GBFF', or 'Unknown'
    """
    try:
        if file_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else: 
            open_func = open
            mode = 'r'
        with open_func(file_path, mode) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    return 'FASTA'
                elif line.startswith('LOCUS'):
                    return 'GBFF'
                else:
                    break  # Only look at the beginning
    except Exception as e:
        print(f"Error reading file: {e}")
        return 'Unknown'

class GBFFSpeciesIdentifier:
    ignored_words = ['Candidatus', 'MAG', 'Uncultured', 'uncultured', 'environmental', 'bacterium']
    pattern = r'((?:' + '|'.join(map(re.escape, ignored_words)) + r')\s+)*([A-Z][a-z]+)\s+([a-z]+)'

    species_regex = re.compile(pattern, re.IGNORECASE)

    def __init__(self, reference_file, history_file):
        self.reference_file = reference_file
        self.history_file = history_file
        self.reference = self._load_reference()
        self.history = self._load_history()

    def _load_reference(self):
        ref = {}
        if os.path.isfile(self.reference_file):
            with open(self.reference_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        ref[row[1]] = row[0]
        return ref

    def _load_history(self):
        hist = {}
        if os.path.isfile(self.history_file):
            with open(self.history_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        hist[row[1]] = row[0]
        return hist

    def _get_aabb(self, genus, species):
        return (genus[:2] + species[:2]).upper()

    def _next_species_code(self, aabb):
        existing = [v for v in self.reference.values() if v.startswith(aabb)]
        nums = [int(v.split('.')[1]) for v in existing if '.' in v]
        next_number = max(nums) + 1 if nums else 1
        return f"{aabb}.{next_number:03d}"

    def _next_strain_code(self, aabb_xxx):
        existing = [v for v in self.history.values() if v.startswith(aabb_xxx + ".")]
        nums = [int(v.split('.')[-1]) for v in existing]
        next_number = max(nums) + 1 if nums else 1
        return f"{next_number:04d}"

    def _extract_info_from_definition(self, record):
        definition = record.description
        print(definition)

        # Step 1: Match genus and species
        match = self.species_regex.search(definition)
        if not match:
            raise ValueError(f"Could not parse genus and species from: {definition}")
        genus = match.group(2)
        species = match.group(3)

        # Step 2: Extract everything after the species
        remainder = definition[match.end():].strip()

        # Step 3: Look for "isolate XYZ" or "strain XYZ"
        strainIDtmp = None
        isolate_match = re.search(r'isolate[_\-\s]+([\w\-\.]+)', remainder)
        strain_match = re.search(r'strain[_\-\s]+([\w\-\.]+(?:[\s_]+[\w\-\.]+)?)', remainder)

        if isolate_match:
            strainIDtmp = isolate_match.group(1)
        elif strain_match:
            strainIDtmp = strain_match.group(1)
        else:
            # fallback: use first word in remainder if any
            strainIDtmp = remainder.split()[0] if remainder else "no_strain"

        # Step 4: Trim bin suffixes like -bin1, _bin.2, .bin3
        strainIDtmp = re.split(r'[-_.]bin[\._]?\d+', strainIDtmp)[0].strip()

        return genus, species, strainIDtmp

    def assign_identifier_from_gbff(self, gbff_path):
        open_func = gzip.open if gbff_path.endswith('.gz') else open
        with open_func(gbff_path, 'rt') as handle:
            record = next(SeqIO.parse(handle, 'genbank'))

        genus, species, strain = self._extract_info_from_definition(record)
        species_key = f"{genus} {species}"
        strain_key = f"{species_key} {strain}"

        # Species-level code
        if species_key not in self.reference:
            aabb = self._get_aabb(genus, species)
            aabb_xxx = self._next_species_code(aabb)
            self.reference[species_key] = aabb_xxx
            with open(self.reference_file, 'a', newline='') as csvfile:
                csv.writer(csvfile, delimiter='\t').writerow([aabb_xxx, species_key])
        else:
            aabb_xxx = self.reference[species_key]

        # Strain-level code
        if strain_key not in self.history:
            yyyy = self._next_strain_code(aabb_xxx)
            identifier = f"{aabb_xxx}.{yyyy}"
            self.history[strain_key] = identifier
            with open(self.history_file, 'a', newline='') as csvfile:
                csv.writer(csvfile, delimiter='\t').writerow([identifier, strain_key])
        else:
            identifier = self.history[strain_key]

        print(f"[INFO] Identifier assigned: {identifier} ({strain_key})")
        return identifier

class SpeciesIdentifier:
    ignored_words = ['Candidatus', 'MAG', 'Uncultured', 'uncultured', 'environmental', 'bacterium']
    pattern = r'((?:' + '|'.join(map(re.escape, ignored_words)) + r')\s+)*([A-Z][a-z]+)\s+(sp\.|[a-z]+)'
    species_regex = re.compile(pattern, re.IGNORECASE)

    def get_identifier_from_file(self, file_path):
        """
        Retrieves only the identifier (AABB.XXX.YYYY) from a FASTA or GBFF file,
        without writing any output files.
        """
        def detect_format(path):
            with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')) as f:
                for line in f:
                    if line.startswith('>'):
                        return 'fasta'
                    elif line.startswith('LOCUS'):
                        return 'gbff'
            raise ValueError("Unknown file format: not FASTA or GBFF.")

        file_format = detect_format(file_path)

        # Extract species info
        if file_format == 'fasta':
            with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
                for line in f:
                    if line.startswith('>'):
                        genus, species, strain = self._extract_species_info(line)
                        break
                else:
                    raise ValueError("No FASTA header found.")
        elif file_format == 'gbff':
            with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
                record = next(SeqIO.parse(f, 'genbank'))
                organism = record.annotations.get("organism", "")
                strain = ""
                for feature in record.features:
                    if feature.type == "source":
                        strain = feature.qualifiers.get("strain", [""])[0]
                        break
                genus, species = organism.split()[:2]

        species_key = f"{genus} {species}"
        strain_key = f"{species_key} {strain}"

        if strain_key in self.history:
            return self.history[strain_key]
        else:
            aabb = self._get_aabb(genus, species)
            aabb_xxx = self.reference.get(species_key)
            if not aabb_xxx:
                aabb_xxx = self._next_species_code(aabb)
                self.reference[species_key] = aabb_xxx
                with open(self.reference_file, 'a', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter='\t')
                    writer.writerow([aabb_xxx, species_key])

            yyyy = self._next_strain_code(aabb_xxx)
            identifier = f"{aabb_xxx}.{yyyy}"
            self.history[strain_key] = identifier
            with open(self.history_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([identifier, strain_key])

            return identifier

    def __init__(self, reference_file, history_file, output_dir):
        self.reference_file = reference_file
        self.history_file = history_file
        self.output_dir = output_dir

        os.makedirs(self.output_dir, exist_ok=True)

        self.reference = self._load_reference()
        self.history = self._load_history()

    def _load_reference(self):
        reference = {}
        if os.path.isfile(self.reference_file):
            with open(self.reference_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        reference[row[1]] = row[0]
        return reference

    def _load_history(self):
        history = {}
        if os.path.isfile(self.history_file):
            with open(self.history_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        history[row[1]] = row[0]
        return history

    def _get_aabb(self, genus, species):
        return (genus[:2] + species[:2]).upper()

    def _next_species_code(self, aabb):
        existing_codes = [
            code for code in self.reference.values()
            if code.startswith(aabb)
        ]
        numbers = [int(code.split('.')[1]) for code in existing_codes if '.' in code]
        next_number = max(numbers) + 1 if numbers else 1
        return f"{aabb}.{next_number:03d}"

    def _next_strain_code(self, aabb_xxx):
        existing = [
            full_id for full_id in self.history.values()
            if full_id.startswith(aabb_xxx + ".")
        ]
        numbers = [int(code.split('.')[-1]) for code in existing]
        next_number = max(numbers) + 1 if numbers else 1
        return f"{next_number:04d}"

    def _extract_species_info(self, fasta_header):
        """
        Extracts genus, species, and remaining description from a FASTA header.
        Uses regex to find the first occurrence of a species name.
        """
        # Remove '>' and sequence ID first
        parts = fasta_header[1:].split(maxsplit=1)
        description = parts[1] if len(parts) > 1 else ""

        match = self.species_regex.search(description)
        if match:
            genus = (match.group(2) or "").strip()
            species = (match.group(3) or "").strip()
            strain = description.replace(match.group(0), '').strip() or "no_strain"
            return genus, species, strain

        raise ValueError(f"Cannot detect species name in header: {fasta_header}")

    def assign_identifier_to_fasta(self, fasta_path):
        # Extract species and strain info from the FIRST header
        if fasta_path.endswith('.gz'):
            with gzip.open(fasta_path, 'rt') as f:
                for line in f:
                    if line.startswith('>'):
                        genus, species, strain = self._extract_species_info(line)
                        break
                    else:
                        raise ValueError("No FASTA header found in file.")
        else :
            with self.open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        genus, species, strain = self._extract_species_info(line)
                        break
                else:
                    raise ValueError("No FASTA header found in file.")
    
        species_key = f"{genus} {species}"            # species-level key
        strain_key = f"{species_key} {strain}"        # full strain key
    
        # Step 1: Get or create species-level code (AABB.XXX)
        if species_key in self.reference:
            aabb_xxx = self.reference[species_key]
        else:
            aabb = self._get_aabb(genus, species)
            aabb_xxx = self._next_species_code(aabb)
            self.reference[species_key] = aabb_xxx
            with open(self.reference_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([aabb_xxx, species_key])
    
        # Step 2: Get or create strain-level code (YYYY)
        if strain_key in self.history:
            identifier = self.history[strain_key]
        else:
            yyyy = self._next_strain_code(aabb_xxx)
            identifier = f"{aabb_xxx}.{yyyy}"
            self.history[strain_key] = identifier
            with open(self.history_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([identifier, strain_key])
    
        # Generate renamed FASTA
        if fasta_path.endswith('.gz'):
            output_fasta = os.path.join(
                self.output_dir,
                os.path.basename(fasta_path).rsplit('.', 2)[0] + '.renamed.fasta'
            )
            with gzip.open(fasta_path, 'rt') as infile, open(output_fasta, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        parts = line[1:].split(maxsplit=1)
                        original_seq_id = parts[0]
                        rest = parts[1] if len(parts) > 1 else ""
                        new_header = f">{identifier} {rest.strip()} |OLD:{original_seq_id}|"
                        outfile.write(new_header + '\n')
                    else:
                        outfile.write(line)
        else:
            output_fasta = os.path.join(
                self.output_dir,
                os.path.basename(fasta_path).rsplit('.', 1)[0] + '.renamed.fasta'
            )
            with open(fasta_path, 'r') as infile, open(output_fasta, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        parts = line[1:].split(maxsplit=1)
                        original_seq_id = parts[0]
                        rest = parts[1] if len(parts) > 1 else ""
                        new_header = f">{identifier} {rest.strip()} |OLD:{original_seq_id}|"
                        outfile.write(new_header + '\n')
                    else:
                        outfile.write(line)
    
        print(f"[INFO] Created: {output_fasta}")
        print(f"[INFO] Identifier for assembly: {identifier} ({strain_key})")
        return identifier, output_fasta
    def assign_identifier_to_file(self, file_path):
        """
        Assigns a unique identifier to a FASTA or GBFF file and creates a renamed FASTA.
        Supports both plain and gzipped files.
        """

        # Determine file format
        def detect_format(path):
            with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')) as f:
                for line in f:
                    if line.startswith('>'):
                        return 'fasta'
                    elif line.startswith('LOCUS'):
                        return 'gbff'
            raise ValueError("Unknown file format: not FASTA or GBFF.")

        file_format = detect_format(file_path)

        # Extract species info
        if file_format == 'fasta':
            with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
                for line in f:
                    if line.startswith('>'):
                        genus, species, strain = self._extract_species_info(line)
                        break
                else:
                    raise ValueError("No FASTA header found.")
        elif file_format == 'gbff':
            with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as f:
                record = next(SeqIO.parse(f, 'genbank'))
                organism = record.annotations.get("organism", "")
                strain = ""
                for feature in record.features:
                    if feature.type == "source":
                        strain = feature.qualifiers.get("strain", [""])[0]
                        break
                genus, species = organism.split()[:2]

        species_key = f"{genus} {species}"
        strain_key = f"{species_key} {strain}"

        # Step 1: Species-level code
        if species_key in self.reference:
            aabb_xxx = self.reference[species_key]
        else:
            aabb = self._get_aabb(genus, species)
            aabb_xxx = self._next_species_code(aabb)
            self.reference[species_key] = aabb_xxx
            with open(self.reference_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([aabb_xxx, species_key])

        # Step 2: Strain-level code
        if strain_key in self.history:
            identifier = self.history[strain_key]
        else:
            yyyy = self._next_strain_code(aabb_xxx)
            identifier = f"{aabb_xxx}.{yyyy}"
            self.history[strain_key] = identifier
            with open(self.history_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([identifier, strain_key])

        # Output file path
        output_fasta = os.path.join(
            self.output_dir,
            os.path.basename(file_path).replace('.gz', '').rsplit('.', 1)[0] + '.renamed.fasta'
        )

        # Create renamed FASTA
        with (gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')) as infile, open(output_fasta, 'w') as outfile:
            if file_format == 'fasta':
                for line in infile:
                    if line.startswith('>'):
                        parts = line[1:].split(maxsplit=1)
                        original_seq_id = parts[0]
                        rest = parts[1] if len(parts) > 1 else ""
                        new_header = f">{identifier} {rest.strip()} |OLD:{original_seq_id}|"
                        outfile.write(new_header + '\n')
                    else:
                        outfile.write(line)
            elif file_format == 'gbff':
                records = SeqIO.parse(infile, "genbank")
                for record in records:
                    original_seq_id = record.id
                    description = record.description
                    header = f">{identifier} {description.strip()} |OLD:{original_seq_id}|"
                    outfile.write(header + '\n')
                    seq = str(record.seq)
                    for i in range(0, len(seq), 60):
                        outfile.write(seq[i:i+60] + '\n')

        print(f"[INFO] Created: {output_fasta}")
        print(f"[INFO] Identifier for assembly: {identifier} ({strain_key})")
        return identifier, output_fasta
    
def run_bakta(fasta_file, output_dir, prefix, threads=4, db_path=None, force=False):
    """
    Run Bakta on the provided FASTA file.
    
    :param fasta_file: Path to input FASTA
    :param output_dir: Directory for Bakta output
    :param prefix: Prefix for output files
    :param threads: Number of threads to use
    :param db_path: Path to Bakta database (if None, uses $BAKTA_DB environment variable)
    """
    if not prefix:
        raise ValueError("Prefix is required for Bakta output naming.")
        
    if db_path is None:
        db_path = os.getenv("BAKTA_DB")
        if db_path is None:
            raise ValueError("No Bakta database provided and BAKTA_DB environment variable not set.")
    if os.path.exists(output_dir):
        if force:
            print(f"[INFO] Removing existing output directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            raise FileExistsError(f"Output directory '{output_dir}' already exists. Use --force to overwrite.")
            
    #os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        "bakta",
        "--db", db_path,
        "--output", output_dir,
        "--prefix", prefix,
        "--threads", str(threads),
        fasta_file 
    ]

    print(f"[INFO] Running Bakta:\n{' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
        print("[INFO] Bakta completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Bakta failed with return code {e.returncode}.")
        raise


def parse_tsv(tsv_file, debug=False):
    features = {}
    file_prefix = None

    with open(tsv_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:  # CHANGED to require at least 9 fields
                continue
            # Fill missing fields with "-"
            while len(parts) < 10:
                parts.append('-')
            seq_id, ftype, start, stop, strand, locus_tag, old_locus_tag, gene, product, dbxrefs = parts[:10]
            old_locus_tag = old_locus_tag.strip()
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

    return features, file_prefix

def generate_new_gene_id(locus_tag, count):
    parts = locus_tag.split(".")
    if len(parts) < 3:
        return f"unknown_{count:05d}"
    aabb = parts[0][:4].lower()
    XXX = parts[1]
    YYYY = parts[2]
    tmp = parts[3][:-1]
    sp = tmp.split("_")
    ZZZZZ = int(sp[1])
    return f"{aabb}{XXX}.{YYYY}.{ZZZZZ}"

   
def build_fasta_header(feat, count_by_type):
    gene = feat['Gene']
    if not gene:
        count_by_type[feat['Type']] += 1
        gene = generate_new_gene_id(feat['SequenceId'], count_by_type[feat['Type']])
    return " ".join([
        feat['LocusTag'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
        gene, feat['OLDLocusTag'], feat['Product'], feat['DbXrefs']
    ])
def build_fasta_header4genes(feat, count_by_type):
    gene = feat['Gene']
    if not gene:
        count_by_type[feat['Type']] += 1
        gene = generate_new_gene_id(feat['SequenceId'], count_by_type[feat['Type']])
    return " ".join([
        feat['SequenceId'], feat['Type'], feat['Strand'], feat['Start'], feat['Stop'],
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

def write_fasta(outfile, header, seq):
    with open(outfile, 'a') as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def process_ffn(input_file, features, file_prefix, suffix_map, out_dirs, count_by_type, debug=False):
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
                        fasta_header = build_fasta_header4genes(feat, count_by_type)
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
                fasta_header = build_fasta_header4genes(feat, count_by_type)
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
                        fasta_header = build_fasta_header4genes(feat, count_by_type)
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
                fasta_header = build_fasta_header4genes(feat, count_by_type)
                write_fasta(outfile, fasta_header, "".join(seq_lines))
                #if debug:
                #    print(f"[DEBUG] Wrote {old_locus} to {outfile}")

def replace_gbff_fields(gbff_file, mapping, locus_tag_map, output_file):
    with open(gbff_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_new_id = None
        in_definition = False
        #key_to_newid, new_id_to_header, old_id_mapping,organism_mapping, replicon_type_mapping
        for line in infile:
            # Replace LOCUS ID with new ID
            if line.startswith("LOCUS"):
                parts = line.split()
                if len(parts) >= 2:
                    old_seqid = parts[1]
                    if old_seqid in mapping['key_to_newid']:
                        current_new_id = mapping['key_to_newid'][old_seqid]
                        new_line = line.replace(old_seqid, current_new_id, 1)
                        outfile.write(new_line)
                        in_definition = False
                        continue
                    else:
                        raise ValueError(f"❌ LOCUS ID '{old_seqid}' not found in FASTA mapping.")

            # Replace DEFINITION with FASTA header
            elif line.startswith("DEFINITION") and current_new_id:
                fasta_header = mapping['new_id_to_header'][current_new_id]
                outfile.write(f"DEFINITION  {fasta_header}\n")
                in_definition = True
                continue

            # Replace ACCESSION with old ID
            elif line.startswith("ACCESSION") and current_new_id:
                old_accession = mapping['old_id_mapping'][current_new_id]
                outfile.write(f"ACCESSION   {old_accession}\n")
                continue

            # Replace VERSION
            elif line.startswith("VERSION") and current_new_id:
                outfile.write("VERSION     Gembases\n")
                continue

            # Replace ORGANISM with extracted one
            elif line.strip().startswith("ORGANISM") and current_new_id:
                organism = mapping['organism_mapping'][current_new_id]
                outfile.write(f"  ORGANISM  {organism}\n")
                continue

            # Replace /mol_type=
            elif '/mol_type=' in line and current_new_id:
                mol_type = mapping['replicon_type_mapping'][current_new_id]
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
            new_locus = parts[0]
            old_locus = parts[6]
            if old_locus and new_locus:
                locus_tag_map[old_locus] = new_locus
    print(f"✅ Built locus_tag mapping for {len(locus_tag_map)} entries from TSV.")
    return locus_tag_map

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
            if row["Type"] == 'oriC':
                continue
            
            if old_seqid not in mapping['key_to_newid']:
                raise ValueError(f"❌ SeqID '{old_seqid}' from TSV not found in FASTA mapping!")

            new_seqid = mapping['key_to_newid'][old_seqid]
            row["#Sequence Id"] = new_seqid

            # Replace strand symbols
            strand = row["Strand"]
            if strand == "+":
                row["Strand"] = "D"
            elif strand == "-":
                row["Strand"] = "C"
            else:
                row["Strand"] = strand

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

    print(f"✅ Temporary TSV written to: {output_file} ({updated_rows} rows updated)")

def parse_args():
    parser = argparse.ArgumentParser(description="Process fasta header and annotation txt file, append taxID using Entrez.")
    parser.add_argument("-f", "--fasta", required=True, help="Input fasta file path")
    parser.add_argument("-t", "--txt", required=True, help="Input annotation txt file path")
    parser.add_argument("-o", "--output", help="Output file path (optional, prints to stdout if omitted)")
    parser.add_argument("-e", "--email", required=True, help="Email address for Entrez API")
    return parser.parse_args()

def parse_fasta_header(header):
    header = header.lstrip('>')
    parts = header.split()
    
    strain_id = parts[0]
    genus = parts[1]
    species = parts[2]
    strain = parts[3] if len(parts) > 3 else ""
    replicons = re.findall(r'\|([^|]+)\|', header)
    
    return strain_id, genus, species, strain, replicons

def get_taxid_from_entrez(genus, species, email):
    Entrez.email=email
    query = f"{genus} {species}"
    handle = Entrez.esearch(db="taxonomy", term=query)
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    if len(id_list) == 0:
        print(f"Warning: No taxID found for {query}", file=sys.stderr)
        return "NA"
    return id_list[0]

def process_txt_file(txt_file_path, strain_id, genus, species, strain_info, replicons, taxid, taxline, output_handle):
    skip_lines = {
        "Bakta:",
        "Software: v1.11.0",
        "Database: v6.0, light",
        "DOI: 10.1099/mgen.0.000685",
        "URL: github.com/oschwengers/bakta"
    }
    
    with open(txt_file_path) as f:
        lines = f.readlines()
    
    for line in lines:
        if line.strip() in skip_lines:
            continue
        if line.startswith("Sequence(s):"):
            replicon_str = " ".join([f"|{r}|" for r in replicons])
            new_line = f"Sequence(s): {strain_id} - {genus} {species} {strain_info} {replicon_str}\n"
            output_handle.write(new_line)
            
            output_handle.write(f"taxID:{taxid}\n")
            output_handle.write(f"taxLine:{taxline}\n")
        else:
            output_handle.write(line)

def get_taxonomy_lineage(taxid,email):
    """
    Given a taxID, fetch taxonomy info from NCBI and return a string:
    "TAXline: Kingdom; Phylum; Class; Order; Family; Genus"
    If some ranks are missing, use 'NA' for those.
    """
    Entrez.email=email
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error fetching taxonomy for taxID {taxid}: {e}")
        return None
    
    if not records or len(records) == 0:
        print(f"No taxonomy record found for taxID {taxid}")
        return None
    
    record = records[0]
    
    # Ranks to extract in order
    desired_ranks = ["kingdom", "phylum", "class", "order", "family", "genus"]
    
    # Initialize dict for rank: name
    rank_dict = {rank: "NA" for rank in desired_ranks}
    
    # The lineage is a list of dicts with 'Rank' and 'ScientificName'
    lineage = record.get("LineageEx", [])
    
    for taxon in lineage:
        rank = taxon.get("Rank", "").lower()
        if rank in desired_ranks:
            rank_dict[rank] = taxon.get("ScientificName", "NA")
    
    # Sometimes the current taxon itself is at one of the desired ranks
    current_rank = record.get("Rank", "").lower()
    current_name = record.get("ScientificName", "NA")
    if current_rank in desired_ranks and rank_dict[current_rank] == "NA":
        rank_dict[current_rank] = current_name
    
    # Format output string
    lineage_str = "; ".join(rank_dict[rank].capitalize() if rank_dict[rank] != "NA" else "NA" for rank in desired_ranks)
    return lineage_str

def process_fasta(file_path, si: SpeciesIdentifier):
    output_path = file_path.rsplit('.', 1)[0] + ".renamed.fasta"

    with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
        current_header = ""
        current_sequence = []

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    # Write previous record
                    outfile.write(current_header + '\n')
                    for seq_line in current_sequence:
                        outfile.write(seq_line + '\n')

                current_sequence = []

                try:
                    parts = line[1:].split(maxsplit=1)
                    original_seqid = parts[0]
                    species_strain = parts[1] if len(parts) > 1 else ""

                    identifier = si.get_identifier_from_fasta_header(line)
                    current_header = f">{identifier} {species_strain} |OLD:{original_seqid}|"
                    print(f"{file_path}: {identifier}")
                except ValueError as e:
                    print(f"Warning in {file_path}: {e}")
                    current_header = line  # Use original if invalid

            else:
                current_sequence.append(line)

        # Write last entry
        if current_header:
            outfile.write(current_header + '\n')
            for seq_line in current_sequence:
                outfile.write(seq_line + '\n')

def check_bakta_output_exists(output_dir, id):
    output_dir = Path(output_dir)
    # Define the *key* output file produced by bakta
    expected_suffixes = [".fna", ".ffn", ".faa", ".tsv", ".txt"]
    expected_files = [output_dir / f"{id}{suffix}" for suffix in expected_suffixes]
    if all(f.exists() for f in expected_files):
        return 1
    else:
        return 0
        

def switch_columns_tsv(input_file, output_file, col1, col2):
    """
    Switches two columns in a TSV file, replaces '#Sequence Id' with 'Replicon' in the header,
    writes the result to a new file, then removes the original input file.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to save the output TSV file.
        col1 (int): Index of the first column to switch (1-based).
        col2 (int): Index of the second column to switch (1-based).
    """
    col1 -= 1  # convert to 0-based indexing
    col2 -= 1

    with open(input_file, 'r', newline='', encoding='utf-8') as infile, \
         open(output_file, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for i, row in enumerate(reader):
            if i == 0:
                # Replace '#Sequence Id' with 'Replicon' in the header
                row = ['Replicon' if cell == '#Sequence Id' else cell for cell in row]
                row = ['#LocusTag' if cell == 'LocusTag' else cell for cell in row]
            if len(row) > max(col1, col2):
                row[col1], row[col2] = row[col2], row[col1]
            writer.writerow(row)

    os.remove(input_file)  # Remove the original file after writing the new one
    print(f"✅ Reformatted TSV written to: {output_file}")

def build_key_mapping(fasta_file):
    """
    Parses a FASTA file to extract key mappings and related metadata.

    Returns:
        dicts:
            key_to_newid           : KEY → new FASTA ID
            new_id_to_header       : new FASTA ID → full FASTA header
            old_id_mapping         : new FASTA ID → OLD ID (from |OLD:XYZ|)
            organism_mapping       : new FASTA ID → organism string
            replicon_type_mapping  : new FASTA ID → replicon type (e.g., 'complete genome')
    """
    key_to_newid = {}
    new_id_to_header = {}
    old_id_mapping = {}
    organism_mapping = {}
    replicon_type_mapping = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description

        # Extract |KEY:...| and |OLD:...|
        key_match = re.search(r"\|KEY:([^\|]+)\|?", header)
        old_match = re.search(r"\|OLD:([^\|]+)\|?", header)

        if not key_match:
            raise ValueError(f"❌ Missing |KEY:XYZ| in header: {header}")

        key = key_match.group(1)
        new_id = header.split()[0]
        key_to_newid[key] = new_id
        new_id_to_header[new_id] = header

        if old_match:
            old_id_mapping[new_id] = old_match.group(1)

        # Parse organism name: use 2nd, 3rd, 4th words of header if available
        header_parts = header.split()
        if len(header_parts) >= 4:
            organism = f"{header_parts[1]} {header_parts[2]} {header_parts[3].rstrip(',')}"
            organism_mapping[new_id] = organism
        else:
            organism_mapping[new_id] = "Unknown"

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
    return {
        "key_to_newid": key_to_newid,
        "new_id_to_header": new_id_to_header,
        "old_id_mapping": old_id_mapping,
        "organism_mapping": organism_mapping,
        "replicon_type_mapping": replicon_type_mapping
    }
def define_log_file(input_file, output_folder):

    #extracts file_id 
    filename = os.path.basename(input_file)
    file_stem = pathlib.Path(filename).stem  # removes .gz or .fna.gz
    file_id = file_stem.split("_")[0]
    #creates log folder
    logs_dir = os.path.join(output_folder, "LOGS")
    os.makedirs(logs_dir, exist_ok=True)
    logfile = os.path.join(logs_dir,f"{file_id}.gembases.log")
    return logfile


def check_gembases_output(logfile,identifier):
    """

    Checks if the pipeline was successful by looking for a specific success message in the log file.



    Parameters:

        log_file (str): Path to the log file.

        entry_name (str): Name of the entry to search for in the success message.

    

    Returns:

        bool: True if the success message is found, False otherwise.
    
    """
    success_phrase1 = f"Gembases entry for '{identifier}' performed successfully"
    success_phrase2 = 'The pipeline was successful'
    try:
        with open(logfile, 'r', encoding='utf-8') as f:
            for line in f:
                if success_phrase1 in line:
                    return True
                elif success_phrase2 in line:
                    return True
            return False
    except FileNotFoundError:
        print(f"❌ Log file not found: {logfile}")
        return False
    
def check_gembases_output_noid(logfile):
    """

    Checks if the pipeline was successful by looking for a specific success message in the log file.



    Parameters:

        log_file (str): Path to the log file.

        entry_name (str): Name of the entry to search for in the success message.

    

    Returns:

        bool: True if the success message is found, False otherwise.
    
    """
    success_phrase1 = "Pipeline already completed Previously!"
    success_phrase2 = 'The pipeline was successful'
    try:
        with open(logfile, 'r', encoding='utf-8') as f:
            for line in f:
                if success_phrase1 in line:
                    return True
                elif success_phrase2 in line:
                    return True
            return False
    except FileNotFoundError:
        print(f"❌ Log file not found: {logfile}")
        return False    
def export_renamed_fasta_from_gbff(gbff_path, identifier, output_fasta):
    """
    Creates a renamed FASTA file from a GBFF(.gz) file using original DEFINITION lines.
    Header format:
    >identifier.suffix original_DEFINITION |KEY:unique_word| |OLD:locus_id|
    """

    # 1. Parse DEFINITIONs
    def parse_definitions(gbff_path):
        open_func = gzip.open if gbff_path.endswith('.gz') else open
        definitions = {}
        with open_func(gbff_path, 'rt') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                definitions[record.id] = record.description
        return definitions

    # 2. Extract KEY from DEFINITION
    def extract_key_from_definition(defn):
        match = re.search(r'(contig\d+|scaffold\d+|chromosome|plasmid\s+\S+)', defn, re.IGNORECASE)
        return match.group(1) if match else "Unknown"

    # Gather definitions
    definitions = parse_definitions(gbff_path)

    suffix_counters = {'c': 1, 'p': 1, 'm': 1, 'd': 1}
    open_func = gzip.open if gbff_path.endswith(".gz") else open

    with open_func(gbff_path, "rt") as handle, open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(handle, "genbank"):
            locus_id = record.id
            defn = definitions.get(locus_id, "")
            key = extract_key_from_definition(defn)
            old_id = locus_id.replace(".1", "")  # Clean trailing version if present

            genome_type = detect_chrom_or_plasmid(defn, defn)
            if genome_type == "chromosome":
                suffix_type = "c"
            elif genome_type == "plasmid":
                suffix_type = "p"
            elif parse_mag_header(defn):
                suffix_type = "m"
            else:
                suffix_type = "d"

            index = suffix_counters[suffix_type]
            suffix = f".{suffix_type}{index:02d}"
            suffix_counters[suffix_type] += 1

            new_id = f"{identifier}{suffix}"
            clean_defn = defn.rstrip(".")

            header = f">{new_id} {clean_defn} |KEY:{key}| |OLD:{old_id}|"
            out_f.write(header + "\n")

            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                out_f.write(seq[i:i+60] + "\n")

    print(f"✅ Renamed FASTA written to: {output_fasta}")

def rename_features_in_gbff(gbff_input, gbff_output, tsv_mapping_file, renamed_fasta):
    """
    Renames:
    - LOCUS name (record.name) from FASTA
    - locus_tag values from TSV

    Leaves:
    - ACCESSION
    - VERSION

    untouched.
    """

    # Step 1: Build OLD → NEW LOCUS name mapping from FASTA
    replicon_mapping = {}
    for record in SeqIO.parse(renamed_fasta, "fasta"):
        header = record.description
        match = re.search(r'\|OLD:([^\|\s]+)', header)
        if not match:
            raise ValueError(f"❌ Missing |OLD:...| in FASTA header: {header}")
        old_id = match.group(1).split('.')[0]       # e.g., NZ_NIGF01000001
        new_locus = header.split()[0]               # e.g., ABUT.001.0002.d01
        replicon_mapping[old_id] = new_locus

    print(f"✅ Loaded {len(replicon_mapping)} replicon mappings from FASTA")

    # Step 2: Load OLD → NEW locus_tag mapping from TSV
    old_to_new_tag = {}
    with open(tsv_mapping_file, 'r') as tsv:
        reader = csv.DictReader(tsv, delimiter='\t')
        for row in reader:
            old = row['OLDLocusTag'].strip()
            new = row['#LocusTag'].strip()
            if old != "-":
                old_to_new_tag[old] = new

    print(f"✅ Loaded {len(old_to_new_tag)} feature tag mappings from TSV")

    # Step 3: Update GBFF
    open_func = gzip.open if gbff_input.endswith(".gz") else open
    with open_func(gbff_input, "rt") as in_handle, open(gbff_output, "w") as out_handle:
        updated_records = []

        for record in SeqIO.parse(in_handle, "genbank"):
            original_accession = record.annotations.get('accessions', [record.id])[0]
            base_id = original_accession.split('.')[0]

            if base_id not in replicon_mapping:
                raise ValueError(f"❌ No renamed LOCUS for: {base_id}")

            new_locus = replicon_mapping[base_id]

            # ✅ Set new LOCUS (name only), but preserve ACCESSION/VERSION
            record.name = new_locus                         # LOCUS line
            record.id = original_accession                  # ACCESSION/VERSION tracking
            if "." in original_accession:
                record.annotations["sequence_version"] = int(original_accession.split(".")[1])  # ✅ restore VERSION
            else:
                record.annotations["sequence_version"] = 1
            # description: reattach using new LOCUS name
            if record.description.startswith(record.id):
                record.description = record.description.replace(record.id, new_locus, 1)
            else:
                record.description = f"{new_locus} {record.description}"

            # ❌ Do NOT change record.annotations['accessions'] or sequence_version

            # ✅ Rename features
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers:
                    old_tag = feature.qualifiers['locus_tag'][0]
                    if old_tag in old_to_new_tag:
                        new_tag = old_to_new_tag[old_tag]
                        feature.qualifiers['locus_tag'][0] = new_tag

                    # ✅ Rename /protein_id and preserve old as /old_protein_id
                    if 'protein_id' in feature.qualifiers:
                        old_prot_id = feature.qualifiers['protein_id'][0]
                        feature.qualifiers['old_protein_id'] = [old_prot_id]
                        feature.qualifiers['protein_id'][0] = new_tag

            updated_records.append(record)

        SeqIO.write(updated_records, out_handle, "genbank")

    print(f"✅ GBFF written: LOCUS updated, ACCESSION/VERSION preserved → {gbff_output}")

def extract_feature_fastas(gbff_path, identifier, output_dir="."):
    """
    Extracts CDS, tRNA, rRNA, and ncRNA features into categorized FASTA files:
    - Genes/{identifier}.nuc
    - Proteins/{identifier}.prt
    - RNA/{identifier}.trna
    - RNA/{identifier}.rrna
    - RNA/{identifier}.ncrna

    Uses feature coordinates for RNA and nucleotide CDS extraction.
    """

    open_func = gzip.open if gbff_path.endswith(".gz") else open

    cds_nuc = []
    cds_prot = []
    trna_seqs = []
    rrna_seqs = []
    ncrna_seqs = []

    # Read GBFF and extract features
    with open_func(gbff_path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seq = record.seq
            for feature in record.features:
                if not feature.location:
                    continue

                ftype = feature.type.lower()
                locus_tag = feature.qualifiers.get("locus_tag", ["unnamed_feature"])[0]

                # CDS: nucleotide + translation (if available)
                if ftype == "cds":
                    nuc_seq = feature.extract(seq)
                    product = feature.qualifiers.get("product", ["-"])[0]
                    header = make_header(locus_tag, "CDS", feature.location, product)
                    cds_nuc.append(SeqRecord(nuc_seq, id=locus_tag, description=header))

                    if "translation" in feature.qualifiers:
                        prot_seq = Seq(feature.qualifiers["translation"][0])
                        cds_prot.append(SeqRecord(prot_seq, id=locus_tag, description=header))

                # tRNA
                elif ftype == "trna":
                    product = feature.qualifiers.get("product", ["-"])[0]
                    header = make_header(locus_tag, "tRNA", feature.location, product)
                    trna_seqs.append(SeqRecord(feature.extract(seq), id=locus_tag, description=header))

                # rRNA
                elif ftype == "rrna":
                    product = feature.qualifiers.get("product", ["-"])[0]
                    header = make_header(locus_tag, "rRNA", feature.location, product)
                    rrna_seqs.append(SeqRecord(feature.extract(seq), id=locus_tag, description=header))

                # ncRNA
                elif ftype == "ncrna":
                    product = feature.qualifiers.get("product", ["-"])[0]
                    header = make_header(locus_tag, "ncRNA", feature.location, product)
                    ncrna_seqs.append(SeqRecord(feature.extract(seq), id=locus_tag, description=header))

    # Create output folders
    os.makedirs(os.path.join(output_dir, "Genes"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "Proteins"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "RNA"), exist_ok=True)

    # Save to FASTA
    def write_records(records, path):
        if records:
            with open(path, "w") as f:
                SeqIO.write(records, f, "fasta")
            print(f"✅ Wrote {len(records)} → {path}")

    write_records(cds_nuc, os.path.join(output_dir, "Genes", f"{identifier}.nuc"))
    write_records(cds_prot, os.path.join(output_dir, "Proteins", f"{identifier}.prt"))
    write_records(trna_seqs, os.path.join(output_dir, "RNA", f"{identifier}.trna"))
    write_records(rrna_seqs, os.path.join(output_dir, "RNA", f"{identifier}.rrna"))
    write_records(ncrna_seqs, os.path.join(output_dir, "RNA", f"{identifier}.ncrna"))

def common_starting_substring(strings):
    if not strings:
        return ""
    shortest = min(strings, key=len)
    for i in range(len(shortest)):
        char = shortest[i]
        for other in strings:
            if other[i] != char:
                return shortest[:i]
    return shortest

def common_starting_words(strings):
    if not strings:
        return ""
    # Split all strings into lists of words
    split_strings = [s.split() for s in strings]
    # Find the shortest list of words
    shortest = min(split_strings, key=len)
    
    common_words = []
    for i in range(len(shortest)):
        word = shortest[i]
        # Check if this word is present at position i in all strings
        if all(words[i] == word for words in split_strings):
            common_words.append(word)
        else:
            break
    
    # Join the common words back into a string
    return " ".join(common_words)

def extract_common_description_prefix_auto(gbff_path):
    descriptions = []
    _, ext = os.path.splitext(gbff_path)
    if ext == '.gz':
        with gzip.open(gbff_path, 'rt') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                descriptions.append(record.description)
    else:
        with open(gbff_path, 'rt') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                descriptions.append(record.description)
    
    return common_starting_words(descriptions)

def summarize_gbff_assembly(gbff_path):
    """
    Parses a GBFF file to extract genome statistics and annotation summary.
    Returns a dictionary with keys:
    - taxID, taxLine, Length, GC, N50, N90, N ratio, coding density, etc.
    """

    open_func = gzip.open if gbff_path.endswith(".gz") else open

    taxid = None
    taxline = None
    total_length = 0
    gc_count = 0
    n_count = 0
    contig_lengths = []
    feature_counts = Counter()
    cds_bp_total = 0
    hypothetical_count = 0

    with open_func(gbff_path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seq = record.seq
            contig_len = len(seq)
            total_length += contig_len
            gc_count += seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
            n_count += seq.upper().count("N")
            contig_lengths.append(contig_len)

            # Taxonomy parsing
            if not taxid:
                for feature in record.features:
                    if feature.type == "source":
                        dbxrefs = feature.qualifiers.get("db_xref", [])
                        for ref in dbxrefs:
                            if ref.startswith("taxon:"):
                                taxid = ref.split(":")[1]
                        taxline = feature.qualifiers.get("organism", [""])[0]
                        break

            # Feature parsing
            for feature in record.features:
                ftype = feature.type.lower()
                feature_counts[ftype] += 1

                if ftype == "cds":
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    cds_bp_total += abs(end - start)

                    # Hypotheticals
                    prod = feature.qualifiers.get("product", [""])[0].lower()
                    if "hypothetical protein" in prod:
                        hypothetical_count += 1

    # Assembly stats
    contig_lengths.sort(reverse=True)
    n50 = compute_nX(contig_lengths, 50)
    n90 = compute_nX(contig_lengths, 90)
    gc_pct = round((gc_count / total_length) * 100, 1) if total_length else 0
    n_pct = round((n_count / total_length) * 100, 1) if total_length else 0
    coding_density = round((cds_bp_total / total_length) * 100, 1) if total_length else 0

    return {
        "taxID": taxid or "unknown",
        "species": taxline or "unknown",
        "taxline": "unknown",
        "Length": total_length,
        "repliconsCount": len(contig_lengths),
        "GC": gc_pct,
        "N50": n50,
        "N90": n90,
        "N ratio": n_pct,
        "coding density": coding_density,
        "tRNAs": feature_counts.get("trna", 0),
        "tmRNAs": feature_counts.get("tmrna", 0),
        "rRNAs": feature_counts.get("rrna", 0),
        "ncRNAs": feature_counts.get("ncrna", 0),
        "ncRNA regions": feature_counts.get("ncrna_region", 0),
        "CRISPR arrays": feature_counts.get("crispr", 0),
        "CDSs": feature_counts.get("cds", 0),
        "pseudogenes": feature_counts.get("pseudogene", 0),
        "hypotheticals": hypothetical_count,
        "sORFs": feature_counts.get("sorf", 0),
        "gaps": feature_counts.get("gap", 0),
        "oriCs": feature_counts.get("origin_of_replication", 0),
        "oriVs": feature_counts.get("oriv", 0),
        "oriTs": feature_counts.get("orit", 0)
    }


def compute_nX(lengths, X):
    """
    Compute N50, N90, etc.
    """
    threshold = sum(lengths) * (X / 100)
    cumulative = 0
    for l in lengths:
        cumulative += l
        if cumulative >= threshold:
            return l
    return 0

def make_header(locus_tag, ftype, loc, product, old_tag=None):
    strand = "D" if loc.strand == 1 else "C"
    start = int(loc.start) + 1
    end = int(loc.end)
    old_part = f" |OLD:{old_tag}|" if old_tag else ""
    return f"{locus_tag} {ftype} {start} {end} {strand} {product}{old_part}"