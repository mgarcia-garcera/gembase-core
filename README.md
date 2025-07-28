# gembase-core
## v.00

this is a repo to store the package gembase-core, a set of scripts orbitting around the "gembase" format, originally ideated by Eduardo Rocha, which I have adopted since my time at Rocha lab. 

I will make evolve this README as I start adding scripts and making them work in this environment.

2025.06.05. I decided that I am converting my scripts to build the gembase + build the core-genome & the pan-genome to python, with the help of MissGPT.

2025.06.06. Connection to git

2025.06.07. Build an awk splitter that takes multi-fastas from RefSeq and splits them into fastas for individual organisms

2025.06.09. Build a python script that takes the header information and reshapes it into a AABB.XXX.YYYY format.

    2025.06.13 - modify the script to improve REGEX of the Species name (it would consider words such as MAG, uncultured, environmental, as part of the species name)

    2025.06.13 - modify the script to solve the issue of the species identifier (XXX) being different for each strain.

2025.06.11. Make bakta work, re-build DB and build a bakta runner py

2025.06.13. Build a bakta Replicon header replacement for a structure that actually works for better parsing.

    2025.06.16 - Modify the script to return the final header in the fasta file

2025.06.13. Build a script that checks GEMBASE directory struture exists and make it if not.

2025.06.16. Build a script that takes the renamed replicon file and the bakta tsv output and returns a re-strucutred tsv with protein identifier following the structure of the gembase identifier.

2025.06.16 Build a script that takes the renamed replicon file and the renamed bakta tsv file and modifies the gbff file, replacing the information where needed.

2025.06.17 Build a script that takes the FFN + FAA and reformats the headers to keep with the AABB.XXX.YYYY format.

2025.06.17 Build a script that takes the replicon fasta and the summary txt and builds the .inf.

2025.06.18 Transfer all functions to GembaseFunctions.py (GF), and build a script called gembases_builder which calls all the functions at GF and smoothly converts a fasta into a gembases.

    2025.06.24 - Modify function parse_tsv to provide Sequence Id as the input for process_ffn and process_faa
    2025.06.24 - Modify function process_ffn and process_faa to create those features where no annotation exists (otherwise, they were skipped)

2025.06.25 Create function in GembaseFunctions.py that parses the LOG file in search of the sentence "Gembases entry for '{AABB.XXX.YYYY}' performed successfully" and stops the pipeline if it finds it.

2025.07.03 Built script that constructs the pangenome for any species identifier "AABB.XXX"
	2025.07.06 - The script works fine when called alone, but on loop it has a bug, which allowed me to see that some species the identifier assigned is not correct.

2025.07.07 Built a script that extracts the information of "Strain identifier, species identifier, species name, TaxID" into a log file
	2025.07.07 - Out of it, I manually built a manually curated list of wrong assignations. (96 species, out of 21193, 0.46%)
2025.07.08 Built a script that corrects a manually curated list of mis-assignations and the correct ones. It outputs the following:

    Processing FAAT.001 -> ZIFA.001 
✅ Found genome IDs: ['FAAT.001.0001']
 Last used index for ZIFA.001 is 0000
FAAT.001.0001 -> ZIFA.001.0001
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Replicons/FAAT.001.0001.fna.gz
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Replicons/ZIFA.001.0001.fna.gz 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Genes/FAAT.001.0001.nuc
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Genes/ZIFA.001.0001.nuc 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Proteins/FAAT.001.0001.prt
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/Proteins/ZIFA.001.0001.prt 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/FAAT.001.0001.trna
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/ZIFA.001.0001.trna 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/FAAT.001.0001.ncrna
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/ZIFA.001.0001.ncrna 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/FAAT.001.0001.rrna
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/RNA/ZIFA.001.0001.rrna 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/FAAT.001.0001.tsv
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/ZIFA.001.0001.tsv 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/FAAT.001.0001.gbff.gz
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/ZIFA.001.0001.gbff.gz 
 updating file: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/FAAT.001.0001.inf
      ✅ renamed to: /scratch/hdd2/mgg/GBFF/Z/Zimmermannella/gembases/LSTINFO/ZIFA.001.0001.inf 


for species with < 200 strains: 
	- Run gembasesBuilderGBFF
	- Run gzipper (to compress GBFF and FNA files
	- Run buildPangenomes4species

For species with >200 strains:
	- skani sketch -t 80 --fast -o skani/$species/$species.skanidb Pre-skani/$capital/$species/*.gz; 
	- skani triangle -t 80 --fast --sparse -o skani/$species/$species.ANI.tb skani/$species/$species.skanidb/*
	- DendogramPlotter.py -i $species.ANI.tb -o $species.clusters.png -c $species.centroids.tsv -d 0.1 -m ward -f png
	- run run_gembases_centroids.py -a $capital/$genus/$species/$species.accessions.lst -c $species.centroids.tsv -r gembases.ref -hi gembases.skani.hist -o $output -e e@mail.com -l $output/LOGS/;
	- run buildPangenomes4species
	 



