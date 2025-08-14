from pathlib import Path
import pathlib
from typing import List
import re
import os
from typing import List
import shutil
import subprocess
import sys
import pandas as pd

def find_genomes_by_taxon(taxon: str, 
                          proteins_dir: Path, 
                          lstinfo_dir: Path) -> List[Path]:
    """
    Find .prt files (proteomes) that belong to a given taxonomic name.
    
    Args:
        taxon (str): Taxonomic name (either gembases format or biological name).
        proteins_dir (Path): Path to the Proteins directory.
        lstinfo_dir (Path): Path to the LSTINFO directory (contains .inf files).
    
    Returns:
        List[Path]: List of matching .prt file paths.
    """
    matched_files = []
    # Case 1: Gembases format (e.g., AABB.XXX)
    if re.fullmatch(r"[A-Z]{4}\.[0-9]{3}", taxon):
        pattern = f"{taxon}.*.prt"
        matched_files = list(proteins_dir.glob(pattern))
    
    else:
        # Case 2: Biological name — search inside .inf files
        taxon = taxon.lower()
        for inf_file in lstinfo_dir.glob("*.inf"):
            with inf_file.open("r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("GEBMASES_Id: "):
                        genome_tmp = line.strip().split()[1]
                    elif line.startswith("taxline:") and taxon in line.lower():
                        genome_id = genome_tmp  # Should be AABB.XXX.YYYY
                        prt_pattern = f"{genome_id}.prt"
                        prt_path = os.path.join(proteins_dir, prt_pattern)
                        matched_files.append(prt_path)
                        break  # Stop after first TAXLINE match
                        
    return matched_files

def count_fasta_headers(filepath: Path) -> int:
    with open(filepath, 'rb') as f:
        return sum(chunk.count(b'>') for chunk in iter(lambda: f.read(8192), b''))

def find_genome_with_most_genes(prt_files: List[Path]) -> str:
    """
    Given a list of .prt files, return the genome_id (filename without extension)
    of the genome with the most protein-coding genes.
    
    Args:
        prt_files (List[Path]): List of .prt files
    
    Returns:
        str: genome_id with the most genes
    """
    max_genes = -1
    best_genome = None

    for prt_file in prt_files:
        gene_count = count_fasta_headers(prt_file)
        genome_id = Path(prt_file).stem
        if gene_count > max_genes:
            max_genes = gene_count
            best_genome = genome_id

    return best_genome

def get_basenames(file_list):
    """
    Convert a list of file paths to a list of basenames (no directory, no extension).

    Args:
        file_list (List[Union[str, Path]]): List of file paths (as str or Path).

    Returns:
        List[str]: List of basenames without extensions.
    """
    basenames = [Path(f).stem for f in file_list]
    return basenames

def parse_opscan_output(stdout):
    results = []
    with open(stdout, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        #awk '/^BB/{escogen=$2; othgen=$4; sim=$6; dif=$7; 
            #getline; othid=$9; escostart=$6; escostop=$7; 
            # getline; 
            # print escogen, othgen, $9,othid, sim, dif, escostart, escostop, $6, $7}' $output/$id1
        if lines[i].startswith("BB"):
            parts1 = lines[i].split()
            escogen, othgen = parts1[1], parts1[3]
            sim, dif = float(parts1[5]), float(parts1[6])

            i+= 1
            parts2 = lines[i].split()
            othdef, escostart, escostop = " ".join(parts2[6:]),parts2[3],parts2[4]

            i+= 1
            parts3 = lines[i].split()
            escodef, othstart, othstop = " ".join(parts3[6:]), parts3[3],parts3[4]
            print(f"{[escogen,othgen,escodef,othdef,sim,dif,escostart, escostop,othstart,othstop]}\n")
            results.append([escogen,othgen,escodef,othdef,sim,dif,escostart, escostop,othstart,othstop])
        i += 1
    
    return results

def save_opscan_log(stderr,logfile):
    with open(logfile, "w") as log:
        for line in stderr:
            log.write(line)

def save_shrt_file(shrt_data, out_path):
    with open(out_path, "w") as f:
        for row in shrt_data:
            f.write("\t".join(map(str, row)) + "\n")

def load_shrt_file(shrt_path):
    shrt_data = []
    print(str(shrt_path))
    with open(str(shrt_path), 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            shrt_data.append(parts)
    return shrt_data

def tsv_format(input_tsv, output_tsv):
    best_hits = {}
    print(f"Parsing {input_tsv}....")
    with open(input_tsv,'r') as input:
        for line in input:
            parts=line.strip().split("\t")
            if len(parts)<4:
                continue

            querydesc=parts[2].split()
            querydesc=querydesc[5:]

            refdesc=parts[3].split()
            refdesc=refdesc[5:]
            parts[2]=" ".join(querydesc)
            parts[3]=" ".join(refdesc)
            query, target, qdesc, tdesc, pident, qcov = parts
            if query not in best_hits:
                best_hits[parts]
            else:
                best_pident = float(best_hits[query][4])
                best_qcov = float(best_hits[query][5])
                if(pident > best_pident) or (pident == best_pident and qcov > best_qcov):
                    best_hits[query] = parts
    print(best_hits)
    with open(output_tsv, 'w') as output:
        for parts in best_hits.values():        
            output.write("\t".join(parts) + "\n")

def log_step(step_desc, result, stdout, stderr):
    header = f"{step_desc}\n\n\n"
    stdout.write(header)
    stdout.write(result.stdout.decode() if isinstance(result.stdout, bytes) else result.stdout)
    stderr.write(header)
    stderr.write(result.stderr.decode() if isinstance(result.stderr, bytes) else result.stderr)

def mmseqs2_pairwise_search(query_file, ref_file, output_tsv, id1, tmpfolder,threads=8):
    workingfolder=os.path.join(tmpfolder,id1)
    os.makedirs(workingfolder)
    try:
        query_db = os.path.join(workingfolder, "queryDB")
        ref_db = os.path.join(workingfolder, "refDB")
        result_db = os.path.join(workingfolder, "resultDB")

        #logs
        log1=os.path.join(tmpfolder,f"{id1}.mmseqs2.log")
        log2=os.path.join(tmpfolder,f"{id1}.mmseqs2.err")
        
        with open(log1,"w") as stdout, open(log2,"w") as stderr:
            print(f"{id1}: CreateDB query")
            r1=subprocess.run(["mmseqs","createdb",query_file,query_db], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check = True)
            log_step(f"STEP1: CreateDB for QUERY GENOME", r1, stdout, stderr)
            print(f"{id1}: CreateDB Reference")
            r2=subprocess.run(["mmseqs","createdb",ref_file,ref_db], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check = True)
            log_step(f"STEP2: CreateDB for REFERENCE GENOME", r2, stdout, stderr)
            print(f"{id1}: search Query VS Reference")
            r3=subprocess.run(["mmseqs","search",query_db, ref_db, result_db, workingfolder, "--threads", str(threads)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check = True)
            log_step(f"STEP3: Search HITS between QUERY AND REFERENCE", r3, stdout, stderr)
            tmp_tsv=os.path.join(tmpfolder,"Tmp.tsv")
            r4=subprocess.run(["mmseqs", "convertalis", query_db, ref_db, result_db, tmp_tsv, "--format-output","query,target,qheader,theader,pident,qcov"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check = True)
            log_step(f"STEP3: print OUTPUT", r4, stdout, stderr)
        
        print(f"Transferring {tmp_tsv} to {output_tsv}")
        tsv_format(tmp_tsv, output_tsv)

    finally:
        shutil.rmtree(workingfolder)
        return(0)

def get_syntenic_stream(lines, memory=5, threshold=10):
    """
    Streaming synteny logic for new tab-delimited format.
    Assumes the 'position' column is the last field (numeric).
    """
    last_positions = []
    for nr, line in enumerate(lines, start=1):
        parts = line.strip().split("\t")
        if not parts:
            continue
        
        try:
            current = float(parts[-1])  # last column as "position"
        except ValueError:
            # If missing numeric position, fallback to row number
            current = nr

        counter = 0
        for prev in last_positions:
            dist = current - prev
            if -threshold <= dist <= threshold:
                counter += 1

        # Maintain rolling memory
        last_positions.append(current)
        if len(last_positions) > memory:
            last_positions.pop(0)

        yield line.strip() + "\t" + str(counter)

def filter_and_synteny_streaming(input_path, output_path, similarity, radius, dist_min, sum_lim, debug=False):
    def debug_print(stage, iterable):
        for item in iterable:
            if debug:
                print(f"[{stage}] {item}")
            yield item

    # Step 1: sort -k 2 (tab-delimited)
    sort_by_col2 = subprocess.Popen(
        ["sort", "-t", "\t", "-k", "2", input_path],
        stdout=subprocess.PIPE,
        text=True
    )

    # Step 2: filter by similarity
    def filter_by_similarity(stream):
        for line in stream:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 5:
                continue
            try:
                sim_val = float(parts[4])
            except ValueError:
                print(f"The value in the identity column is not numeric: {parts[4]}")
            if sim_val >= similarity:
                yield "\t".join(parts)

    filtered_gen = debug_print("filter_by_similarity", filter_by_similarity(sort_by_col2.stdout))

    # Step 3: add NR
    def add_nr(rows):
        for nr, row in enumerate(rows, start=1):
            yield row + "\t" + str(nr)

    numbered_gen = debug_print("add_nr", add_nr(filtered_gen))

    # Step 4: sort -k 1
    sort_by_col1 = subprocess.Popen(
        ["sort", "-t", "\t", "-k", "1"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True
    )
    for row in numbered_gen:
        sort_by_col1.stdin.write(row + "\n")
    sort_by_col1.stdin.close()

    # Step 5: get_syntenic
    synt_gen = debug_print(
        "get_syntenic_stream",
        get_syntenic_stream(sort_by_col1.stdout, memory=radius, threshold=dist_min)
    )

    # Step 6: filter by lim
    def filter_lim(rows):
        for line in rows:
            parts = line.strip().split("\t")
            if int(parts[-1]) >= sum_lim:
                yield "\t".join(parts)

    final_gen = debug_print("filter_lim", filter_lim(synt_gen))

    # Step 7: sort -k 1 again
    sort_final = subprocess.Popen(
        ["sort", "-t", "\t", "-k", "1"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True
    )
    for row in final_gen:
        sort_final.stdin.write(row + "\n")
    sort_final.stdin.close()

    # Step 8: write output
    with open(output_path, "w") as out:
        for line in sort_final.stdout:
            if debug:
                print(f"[final_sort] {line.strip()}")
            out.write(line)
            
def run_mmseqs_and_parse(listfiles, reference, input_folder, output_folder, tmpfolder, sum_lim, radius, dist_min, synteny, identity):

    #define things that are not defined
    reference_file = str(f"{input_folder}/{reference}.prt")
    for file2 in listfiles:
        element2 = Path(file2).stem
        if element2 == reference:
            continue

        id1 = f"{reference}.{element2}.opsc"
        id2 = f"{reference}.{element2}.opsc.{synteny}.{identity}"
    
        print(f"[Processing pair] {reference} vs {element2}")
        shrt_file_tmp = Path(f"{tmpfolder}/{id1}.shrt")
        shrt_file_out = Path(f"{output_folder}/{id1}.shrt")

        if not shrt_file_tmp.exists():
            print(f"     STEP {reference} vs {element2}, not yet performed") 
            mmseqs2_pairwise_search(file2,reference_file,shrt_file_out,id1,tmpfolder,threads=40)
            shutil.copy(shrt_file_out, tmpfolder)    
        else:
            print(f"     STEP {reference} vs {element2}, already performed")
            shutil.copy(shrt_file_tmp, shrt_file_out)

        # ----- step 2: synteny ------

        synt_file_out = Path(f"{output_folder}/{id2}.synt")
        synt_file_tmp = Path(f"{tmpfolder}/{id2}.synt")

        if not synt_file_out.exists():
            filter_and_synteny_streaming(shrt_file_out, synt_file_out, identity, radius, dist_min, sum_lim)
            shutil.copy(synt_file_out, tmpfolder)
        else:
            print(f"     STEP Synteny {reference} vs {element2}, already performed")
            shutil.copy(synt_file_tmp, synt_file_out)
    return(0)
    
def filter_and_synteny(shrt_data, similarity, radius, dist_min, sum_lim):
    """
    Replacement for sort/awk pipeline to detect syntenic blocks
    This is a simplified logic mirroring the filtering and scoring steps
    """
    #1. Filter by similarity
    filtered = [row for row in shrt_data if row[4]>=similarity]

    #2. sort by second column
    filtered.sort(key=lambda r: r[1])

    #3. add running index (like awk's NR)
    for idx, row in enumerate(filtered, start = 1):
        row.append(idx)

    #4. sort by first column
    filtered.sort(key=lambda r: r[0])

    #5. identify synteny blocks
    synt_blocks = []
    for idx, row in enumerate(filtered) :
        count_nearby = sum(
            1 for other in filtered
            if abs(int(row[-1]) - int(other[-1])<= radius)
        )
        if count_nearby >= dist_min:
            score = count_nearby + int(row[-1])
            if score >= sum_lim:
                synt_blocks.append(row[:6] + [score])

    synt_blocks.sort(key=lambda r: r[0])
    return synt_blocks

def save_synt_file(synt_data, out_path):
    with open(out_path, "w") as f:
        for row in synt_data:
            f.write("\t".join(map(str,row))+"\n")

def run_opscan_and_parse(listfiles, reference, mm, input_folder, output_folder, tmpfolder,sum_lim, radius,dist_min,synteny,identity):

    #define things that are not defined
    reference_file = str(f"{input_folder}/{reference}.prt")
    for file2 in listfiles:
        element2 = Path(file2).stem
        if element2 == reference:
            continue

        id1 = f"{reference}.{element2}.opsc"
        id2 = f"{reference}.{element2}.opsc.{synteny}.{identity}"
    
        print(f"[Processing pair] {reference} vs {element2}")
        shrt_file_tmp = Path(f"{tmpfolder}/{id1}.shrt")
        shrt_file_out = Path(f"{output_folder}/{id1}.shrt")

        if not shrt_file_tmp.exists():
            print(f"     STEP {reference} vs {element2}, not yet performed") 
            cmd = [
                  "opscan", "-H", 
                  "-M", str(mm), 
                  "-t","37", 
                  "-r", "1.3", 
                  "-F", "-E", 
                  "-c", "-Q", 
                  "-U", "-O", 
                  str(reference_file), str(file2)
             ]
            #with open (shrt_file_out, "w") as out_f, open (f"{shrt_file_out}.log", "w") as log_f:
            #    subprocess.run(cmd, stdout = out_f, stderr = log_f, check = True)

            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            with open(f"{shrt_file_out}", "w") as f:
                f.write(result.stdout)
            with open(f"{shrt_file_out}.log", "w") as f:
                f.write(result.stderr)
            #check if the run
            print("Return code", result.returncode)

            #store the data modified in file
            shrt_data = parse_opscan_output(shrt_file_out)
            save_shrt_file(shrt_data,shrt_file_out)


            print("STDERR:")
            save_opscan_log(result.stderr, f"{tmpfolder}/{id1}.log")
            shutil.copy(shrt_file_out, tmpfolder)    
        else:
            print(f"     STEP {reference} vs {element2}, already performed")
            shutil.copy(shrt_file_tmp, shrt_file_out)

        # ----- step 2: synteny ------

        synt_file_out = Path(f"{output_folder}/{id2}.synt")
        synt_file_tmp = Path(f"{tmpfolder}/{id2}.synt")

        if not synt_file_out.exists():
             shrt_data = load_shrt_file(shrt_file_out)
             synt_data = filter_and_synteny(shrt_data, identity, radius, dist_min, sum_lim)

             save_synt_file(synt_data, synt_file_out)

             
             shutil.copy(synt_file_out, tmpfolder)
        else:
            print(f"     STEP Synteny {reference} vs {element2}, already performed")
            shutil.copy(synt_file_tmp, synt_file_out)
    return(0)
    
def read_table(filepath):
    """Read a whitespace-delimited file into a list of rows (list of lists)."""
    with open(filepath) as f:
        return [line.strip().split() for line in f if line.strip()]
    
def write_table(filepath, rows):
    """Write rows (list of lists) to a whitespace-delimited file."""
    with open(filepath, "w") as f:
        for row in rows:
            f.write(" ".join(map(str, row)) + "\n")

def join_tables(table1, table2, key_col1=0, key_col2=0):
    """
    Perform an inner join on two tables based on the key column.
    key_col1 and key_col2 are zero-based column indices.
    Returns a new table with merged columns.
    """
    # Build dictionary from table2
    dict2 = {row[key_col2]: row for row in table2}

    joined = []
    for row1 in table1:
        key = row1[key_col1]
        if key in dict2:
            joined.append(row1 + dict2[key])
    return joined

def orthoFinder(workingfolder, taxon,identity,synteny):
    synt_files = list(Path(workingfolder).rglob("*.synt"))
    filenames = [ f.name for f in synt_files ]
    n = len(filenames)
    print(f"There are {n} files with the selected name")
    if n == 0:
         print("No .synt files found in folder")
         sys.exit()
    if n == 1:
        table = read_table(Path(f"{workingfolder}/{filenames[0]}"))
        extracted = [[row[0], row[2], row[1], row[4]] for row in table]
        write_table(f"{workingfolder}/CoreGenome-{taxon}.lst", extracted)
        return(0)
    else:
        # Join first two files
        f1 = os.path.join(workingfolder, f"{filenames[0]}")
        f2 = os.path.join(workingfolder, f"{filenames[1]}")
        t1 = pd.read_csv(f1, sep="\t", header=None)
        t2 = pd.read_csv(f2, sep="\t", header=None)
        
        joined = pd.merge(t1,t2, on=1)  # default join on col0 of both
        joined.columns = range(joined.shape[1])
        joined = joined[[1,2,0,4,8,11]]
        joined.columns = range(joined.shape[1])
        joined = joined.loc[joined.groupby(0)[joined.columns[-1]].idxmax()]
        nfiles=2
        print(f"Core-genome construction for {taxon}: {nfiles} files added: {len(joined)} genes in the Core-genome")
        # Iteratively join with remaining files
        for i in range(2, n):
            nfiles += 1
            f_next = os.path.join(workingfolder, f"{filenames[i]}")
            t_next = pd.read_csv(f_next, sep="\t",header=None)
            t_next = t_next[[1,0,4]]
            t_next.columns = range(t_next.shape[1])
            joined = pd.merge(joined,t_next,on=0)
            joined.columns = range(joined.shape[1])
            joined = joined.dropna()
            joined = joined.loc[joined.groupby(0)[joined.columns[-1]].idxmax()]
            joined.columns = range(joined.shape[1])
            print(f"Core-genome construction for {taxon}: {nfiles} files added: {len(joined)} genes in the Core-genome")
            
            
        joined.to_csv(os.path.join(workingfolder, f"CoreGenome-{taxon}.{identity}.{synteny}.lst"),index=False, header=False, sep="\t")        
        return(0)
