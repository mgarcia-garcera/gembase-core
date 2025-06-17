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


