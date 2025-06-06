#!usr/bin/perl

#this script generates a GEMBASES type of folder, with all the specific data, associated to a genome or a draft.
#It is meant to follow my pipeline for genome assembly and annotation.

use strict;
#get all four variables

chomp( my $folder =$ARGV[0]);
chomp( my $species =$ARGV[1]);
if ($species eq 'Multiple_species') {
    $species = 'Test_stet';
}else{
}

my $genus = $1 if $species=~/^([^\_]+)/;
chomp( my $out_folder =$ARGV[2]);
mkdir $out_folder unless -d $out_folder;
chomp( my $force = $ARGV[3]);

if ($force == 1) {
    qx/rm -r $out_folder/ if -d $out_folder;
}
#undefine $folder if the Input folder does not exist.
if ( ! -d $folder ) {
    undef $folder;
}


if (!$folder || ! $species || ! $out_folder || $folder eq '-h' || $folder eq '--h' || $folder eq '-help' || $folder eq '--help') {
    print "usage: perl $0 <Input Folder> <Species> <Output Folder> <Force>\n";
    print " -h this message\n";
}

#Check if the $folder/PROKKA does exist,  and if the $folder/spades does exist.
if(!-d "$folder/PROKKA" ){
    die "This genome hasn't been yet assembled or the format is not correct. Aborting\n";
}
if(!-d "$folder/spades"){
    die "This genome hasn't been yet assembled or the format is not correct. Aborting\n";
}

my %taxnames;
#construct the out_folder;
die "out_folder exists" if -d $out_folder;
mkdir $out_folder;
mkdir $out_folder . "/LSTINFO/";
mkdir $out_folder . "/Genes/";
mkdir $out_folder . "/Replicons/";
mkdir $out_folder . "/Proteins/";
mkdir $out_folder . "/RNA/";



###extract script location;
chomp(my $location = $0);
my @loc = split /\//, $location;
pop @loc;
$location = join "/", @loc;
if (!$location) {
    chomp($location = qx/pwd/);
}
$location .= '/';
my $gc = $location . 'gc.prt';
my $nodes = $location . 'nodes.bin';
my $names = $location . 'names.list';
our %GC = gencode($gc);
my $species_id = lc($1) . $2 if $species=~/^(\S{2})\S+[\s\_](\S{2})\S+/;
#define the important files:

opendir REP, $folder . 'PROKKA/';
my @rep = grep /\.fsa/, readdir REP;
my  $rep = $folder . 'PROKKA/' . $rep[0];
opendir GEN, $folder . 'PROKKA/';
my @tbl = grep /\.tbl/, readdir TBL;
my  $tbl = $folder . 'PROKKA/' . $tbl[0];
opendir GFF, $folder . 'PROKKA/';
my @gff = grep /\.gff/, readdir GFF;
my  $gff = $folder . 'PROKKA/' . $gff[0];



closedir REP;
closedir GEN;
closedir PRT;
closedir  GFF;
closedir TBL;
undef @tbl;
undef @rep;
undef @gff;
#define the file identifier
my $i = 1;
my $zero1 = 3;
my $name = $species_id;
my $remaining = $zero1 - length($i);
$name .= "." . "0"x$remaining . $i. ".d01.00";
my $identifier = uc($species_id) . "0"x$remaining . $i . "d01";
my %genome_data;
#Start with replicon data
open(GNM, $rep);
open(OUT, ">", $out_folder . "/Replicons/" . $name . ".fst");
    my $header;
    my $zero2 = 5;
    my $j = 1;
    my %contigs;
    my $size;
    while (<GNM>) {
        if (/>(\S+)/) {
            $header = $1;
            my $dif = $zero2 - length($j);
            my $new_header = ">" .  $identifier . "_CONTIG_" . "0"x$dif . $j . " " . $header . "\n";
            $contigs{$header}  = $j;
            print OUT $new_header;
            $j++;
        }else{
            print OUT $_;
            chomp;
            $genome_data{$header} .= $_;
            $size += length($_);
        }
        
    }
    close GNM;
    close OUT;
#now extract the data from tbl;
my $id = 'DRAFT';
        my $prt = $out_folder . '/Proteins/' . $name . ".prt";
        my $gen = $out_folder . '/Genes/' . $name . ".gen";
        my $lst = $out_folder . '/LSTINFO/' . $name . ".lst";
        my $trna = $out_folder . '/RNA/' . $name . ".trna";
        my $rrna = $out_folder . '/RNA/' . $name . ".rrna";
        my $mrna = $out_folder . '/RNA/' . $name . ".tmrna";
        open(TABLE, $tbl);
        open(LST, ">", $lst);
        open(PRT, ">", $prt);
        open(GEN, ">", $gen);
        open(TRNA, ">", $trna);
        open(RRNA, ">", $rrna);
        open(MRNA, ">", $mrna);
        my ($contig,$beg,$end,$product,$contigID);
        my ($EC,$gene,$tag,$prod,$uniprot);
        my $gene_ID;
        my $k = 1;
        my $zero3 = 5;
        while (<TABLE>) {
            if (/>Feature\s(\S+)/) {
                $contig = $1;
                my $tmp= $contigs{$contig};
                my $dif = $zero2 - length($tmp);
                $contigID = 'CONTIG_' . "0"x$dif . $tmp;
            }elsif(/^(\d+)\s+(\d+)\s+(\S+)/){
                undef $EC;
                undef $gene;
                undef $tag;
                undef $prod;
                ($beg,$end,$product) = ($1,$2,$3);
                my $size = abs($end - $beg) + 1;
                my $dif = $zero3 -length($k);
                $gene_ID = $identifier . "_" . "0"x$dif . $k . 0;
                $k++;
            }else{
                if (/EC_number\s+(\S+)/) {
                    $EC = $1;
                }elsif(/gene\s+(\S+)/){
                    $gene = $1;
                }elsif(/locus_tag\s+(\S+)/){
                    $tag = $1;
                }elsif(/UniProtKB:(\S+)/){
                    $uniprot = $1;
                }elsif(/product\s+(.*)\n/){
                    $prod = $1;
                    my $seq;
                    my $direction;
                    if ($beg > $end) {
                        ($beg,$end) = ($end,$beg);
                        $seq = substr($genome_data{$contig},($beg-1),($end+1-$beg));
                        $seq = reverse($seq);
                        $seq =~tr/ATGC/TACG/;
                        $direction  = 'C';
                    }else{
                        $seq = substr($genome_data{$contig},($beg-1),($end+1-$beg));
                        $direction = 'D';
                    }
                    my @seq = unpack "(A3)*", $seq;
                    my $first = $seq[0];
                    my $last = $seq[scalar(@seq) - 1];
                    my $LST = $contigID . " " . $beg . " " . $end . " " . $direction . " " . $product . " " . $gene_ID . " Valid " ;
                    my $GN = $gene_ID . " " . $direction . " " . $first . " " . $last . " " . $beg . " " . $end . " Valid "; 
                    if ($gene) {
                        $LST .= $gene;
                        $GN .= $gene; 
                    }else{
                        $tag = $1 if $tag=~/\_(\d+)/;
                        my $head = lc($1) if $id=~/^(\S{7})/;
                        $tag = $head . "_" . $tag;
                        $LST .= $tag;
                        $GN .= $tag;
                    }
                    $LST .= " | " . $contig;
                    $GN .= " | " . $contig ." ". $contigID;
                    if ($uniprot) {
                        $LST .= " | " . $uniprot;
                        $GN .= " | " . $uniprot;
                    }else{
                        
                    }
                    
                    $LST .= " | " . $prod;
                    $GN .= " | " . $prod;
                    print LST $LST . "\n";
                    $seq = join "\n", unpack "(A60)*", $seq;
                    
                    if ($product eq 'CDS') {
                        $_ = &codon($_,11) for @seq;
                        pop @seq;
                        my $prot = join "", @seq;
                        $prot = join "\n", unpack "(A60)*", $prot;
                        print PRT ">" . $GN . "\n" . $prot . "\n";
                        print GEN ">" . $GN . "\n" . $seq . "\n";
                    }elsif($product eq 'tRNA'){
                        print TRNA ">" . $GN . "\n" . $seq . "\n";
                    }elsif($product eq 'rRNA'){
                        print RRNA ">" . $GN . "\n" . $seq . "\n";
                    }elsif($product eq 'tmRNA'){
                        print MRNA ">" . $GN . "\n" . $seq . "\n";
                    }else{
                        my $wow;
                    }
                }
                
            }
            
        }
        close TABLE;
        close LST;
        close GEN;
        close PRT;
        close RRNA;
        close TRNA;


exit;


sub codon{
    my ($triplet, $code) = @_;
    my $aa = $GC {$code}{code}{$triplet};
    return $aa;
}
sub br{
    my ($gi,$file)=@_;
    my ($pos,$counter,$pos2,$bin_val, $real_val)=0;
    my $buffer;
    open(TAX_ID, $file) || die "$0";
    $pos2=$gi*4;
    binmode TAX_ID;
        
    seek TAX_ID, $pos2,0;
    read(TAX_ID,$buffer,4);
    $real_val=unpack("N",$buffer);
    close TAX_ID;
    return $real_val;
}
sub gencode{
    my ($file) = @_;
    open F, $file;
    my %out;
    while(<F>){
        next if /^\-\-/;
        if (/^\s+name\s\"(.*)\"\s,/) {
            my $name = $1;
            my $line = <F>;
            if ($line=~/^\s+name\s\"(.*)\"\s,/) {
                $line =<F>;
            }
            my $id = $1 if $line =~/^\s+id\s(\d+)\s,/;
            do{
                $out {$id}{name} = $name;
                $line = <F>;
                my $trans = $1 if $line=~/\s+ncbieaa\s+\"(.*)\"\,/;
                $line = <F>;
                my $start = $1 if $line =~/sncbieaa "(.*)"/;
                $line = <F>;
                my $base1 = $1 if $line =~/\s+\-\-\sBase1\s+([ATCG]+)\n/;
                $line = <F>;
                my $base2 = $1 if $line =~/\s+\-\-\sBase2\s+([ATCG]+)\n/;
                $line = <F>;
                my $base3 = $1 if $line =~/\s+\-\-\sBase3\s+([ATCG]+)\n/;
                my @trans = split //, $trans;
                my @start = split //, $start;
                my @base1 = split //, $base1;
                my @base2 = split //, $base2;
                my @base3 = split //, $base3;
                while (@trans) {
                    my $aa = shift @trans;
                    my $start = shift @start;
                    my $triplet = shift @base1;
                    $triplet .= shift @base2;
                    $triplet .= shift @base3;
                    if ($aa eq '*') {
                        $out {$id} {code} {$triplet} = '*';
                        push @{$out {$id}{trans}{STOP}}, $triplet;
                    }else{
                        push @{$out {$id}{trans}{$aa}}, $triplet;
                        $out {$id} {code} {$triplet} = $aa;
                    }
                    if ($start eq '-') {
                        $out {$id}{start}{$triplet} = 0;
                    }else{
                        $out {$id}{start}{$triplet} = 1;
                    }
                }
            }unless !$id;
        }
        my $wow = 'Wow';
    }
    return %out;

}