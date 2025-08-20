#Note: this is pseudocode with the correct command lines used. 
#But the directions are not correct, the same way as the 
#locations of scripts are not correct. This is just for me 
#to keep track of how the method works.

#step 0: Define variables:
assembly_file=assembly_summary.l200.final.tsv
GBFF_folder=output.folder
LOG1=NCBIDump.download.log
hist=gembases.hist
ref=gembases.ref
#step 0: download RefSeq bacteria
    #NCBIRefSeq_DUMP.py uses the assembly_summary.tsv which 
    #contains a column "Capital" and a column "Genus", as well 
    #as the correct location of the genome in the ftp, and 
    #downloads it in batches of XXX genomes at a time.
python NCBIRefSeq_DUMP.py -assembly $assembly_file -o $GBFF_folder -l $LOG1 -s gbff --batchsize 5000
#step 1: build gembasesDB
            #              capital letter of the genus
            #                Genus 
            #                  Species 
for folder in $GBFF_folder/*/*/*; do # this could be parallelized per capital
    #define genus and capital by extracting them from path
    species=$(basename $folder); 
    lst=$(ls $folder | wc -l); 
    $genus=$(dirname $folder)
    capital=$(dirname $genus)
    $GBFF_main=$(dirname $capital)
    genus=$(basename $genus)
    capital=$(basename $capital)

    if [[ "$lst" -gt 100 ]]; then 
        if [[ "$lst" -lt 5000 ]]; then #This is because >5000 was interfering with usual GEN work. Otherwise it could have been done in the same pipeline
            echo $species"    "$lst; 
            mkdir $folder/skani; 
            skani sketch -t 20 --fast -o $folder/skani/$species.skani.db $folder/*.gz; 
            skani triangle -t 20 --fast --sparse -o $folder/skani/$species.ANI.tbl $folder/skani/$species.skani.db/*;
            #dendroplotter is a python script that builds a dendrogram, and extracts the centroids out of skani with 99% identity and 99% coverage
            python DendogramPlotter.py -i $folder/skani/$species.ANI.tbl -o $folder/skani/$species.plot.png -c $folder/skani/$species.centroids.tsv -m ward -f png
            #run_gembases_centroids.py extracts the identifiers from the centroids and runs gembases_builderGBFF.py out of them. It uses a file ".accessions.lst, which only links accessions to files, it is just built with grep"
            python run_gembases_centroids.py -a $folder/skani/$species.accessions.lst \
                                             -c $folder/skani/$species.centroids.tsv \
                                             -r gembases.ref -hi gembases.hist.sort \
                                             -o $GBFF_main/$capital/$genus/gembases/ \
                                             -l $GBFF_main/$capital/$genus/gembases/LOGS/ -e mgarcia@toto.org

        else
            echo $species" Will be done in Obelix"
        fi; 
    else:
        for file in $folder/*.gbff; do
            #define the filename. This will only be used for the log
            filename=$(basename $filename)
            filename=${filename%genomic*}
            #build the output folder. It could be implemented internally but for now it is external.
            if [ ! -d $GBFF_main/$capital/$genus/gembases/ ]; then 
                mkdir $GBFF_main/$capital/$genus/gembases/
            fi
            if [ ! -d $GBFF_main/$capital/$genus/gembases/LOGS/ ]; then
                mkdir $GBFF_main/$capital/$genus/gembases/LOGS/
            fi

            python gembases_builderGBFF.py -i $file -r $ref -hi $hist -o $GBFF_main/$capital/$genus/gembases/ \
                                           -t 40 -e mgarcia@toto.org -l $GBFF_main/$capital/$genus/gembases/$filename.gembases.log
    fi; 
done
#the gembasesDB structure is contrained inside the gembases folder and contains 5 subfolders and 10 files.
#all files belonging to one original GBFF have the same prefix: AABB.XXX.YYYY.{suffix} where:
    #AABB is the first 2 letters from Genus(AA) and first 2 letters from species (BB).
    #XXX is a disambiguation identifier, for when Genome A and genome B have the same combination 
        #of AA BB but their species is different. i.e. ABPR.001: Abyssibacter profundi and ABPR.002:Candidatus Absconditicoccus praedator
    #YYYY is a genome identifier. In theory, the script allows to go up to 99999 genomes. In practice, rarely we go further than 10000 genomes.
#the folder structure and file suffixes are as follows:
    #gembases/Replicons/ which contains .fna.gz files (replicon files)
    #gembases/Genes/ which contains .nuc files (gene multi nucleotide fasta files)
    #gembases/Proteins/ which contains .prt files (CDS multifasta protein files)
    #gembases/RNA/ which contains .rrna, .trna, .ncrna files (functional RNA files in nucleotide multifasta)
    #gembases/LSTINFO/ which contains: .gbff.gz -> GBK file reannotated; 
                                     # .inf     -> contains general information about the genome (taxID, name, tax line, assembly stats, etc.)
                                     # .tsv     -> genome annotation in Tab separated format
#step 2: build Pangenome
for file in $GBFF_folder/*/*/gembases/LSTINFO/*0001.inf; #the filename by choosing the YYYY = 0001 we are choosing the first genome in the species, and hence we get all of them.
    do speciesID=$(basename $file | cut -d'.' -f1,2); 
        capital=${speciesID:0:1}; 
        echo $speciesID"   "$capital;  
        folder=$(dirname $file); 
        folder=$(dirname $folder); 
        #building coding PanDB
        python /scratch/hdd3/mgg/notebooks/labbooks/Marc/Gembases2NextflowGembases/gembase-core/BuildPangenome4species.py -s $speciesID \
            -f $folder -m mmseqs -o /scratch/hdd3/mgg/Refseq_bacteria/panbases/$capital/ \
            -i 0.95 -c 0.95 -l /scratch/hdd3/mgg/Refseq_bacteria/panbases/LOGS/$speciesID.log \
            --tmp /scratch/hdd3/mgg/Refseq_bacteria/tmp/tmp$capital/ --tax /scratch/hdd2/mgg/ListOfTaxids ; 
        #building nonCoding PanDB
        folder=$folder/RNA/
        python /scratch/hdd3/mgg/notebooks/labbooks/Marc/Gembases2NextflowGembases/gembase-core/BuildPangenome4species.py -s $speciesID \
                        -f $folder -m mmseqs -o /scratch/hdd3/mgg/Refseq_bacteria/panbases_noncoding/$capital/ \
                        -i 0.95 -c 0.95 -l /scratch/hdd3/mgg/Refseq_bacteria/panbases_noncoding/LOGS/$speciesID.log \
                        --tmp /scratch/hdd3/mgg/Refseq_bacteria/tmp/tmp$capital/ --tax /scratch/hdd2/mgg/ListOfTaxids ; 
       done
#step 3: concatenate
for capital in {A..Z}; do
    cat /scratch/hdd3/mgg/Refseq_bacteria/panbases/$capital/*.fna > /data/databases/pandb/$capital.fna
    cat /scratch/hdd3/mgg/Refseq_bacteria/panbases_noncoding/$capital/*.fna > /data/databases/pandb/$capital.nc.fna
done
