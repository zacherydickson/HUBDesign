#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Class::Struct;
use Bio::SeqIO;

sub rc($);

struct (GENE => {
            name => '$',
            txid => '$',
            id => '$',
            version => '$',
            chr => '$',
            start => '$',
            end => '$',
            strand => '$'
    });

if(@ARGV < 3) {
    die "Usage: ".basename($0)." OutputDirectory GeneInfo InputFile1 ...\n";
}

my $OutDir = shift(@ARGV);
my $GeneInfoFile = shift(@ARGV);

##Load Gene Info
my %GenesInSeq; #Hash with keys of seqid and values of array ref with elements of GENE Obj
#my %OutputObjects; #Hash with keys of gene name and values of Bio::SeqIO Objects
if(open(my $fh,$GeneInfoFile)){
    while(my $line = <$fh>){
        chomp($line);
        my ($name, $txid, $id, $ver, $chr, $start, $end, $strand) = split(/\t/,$line);
        $GenesInSeq{$chr} = [] unless(exists $GenesInSeq{$chr});
        push(@{$GenesInSeq{$chr}},GENE->new(name => $name, txid => $txid, id => $id, version => $ver, chr => $chr, start => $start, end => $end, strand => $strand)); 
    }
    close($fh);
} else {
    die "Could not open $GeneInfoFile: $!\n"
}

##Process SequenceFiles
while(my $file = shift(@ARGV)){
    if($file =~ m/\.gz$/){
        $file = "zcat $file |";
    }
    open (my $fh, $file) or die "Could not open $file: $!\n";
    my $IN = new Bio::SeqIO(-fh => $fh, -format => 'fasta');
    while(my $seqObj = $IN->next_seq){
        my $acc = $seqObj->id;
        if(exists $GenesInSeq{$acc}){
            #Send each gene to its appropriate file
            foreach my $GeneObj (@{$GenesInSeq{$acc}}){
                my $seq = $seqObj->subseq($GeneObj->start,$GeneObj->end);
                $seq = rc($seq) if($GeneObj->strand == -1);
                my $outseqObj = Bio::Seq->new(-seq => $seq, -id => $GeneObj->id);
                my $OUT = Bio::SeqIO->new(-file => ">>$OutDir/".$GeneObj->name.".fna",-format => 'fasta');
                $OUT->write_seq($outseqObj);
            }
        }
    }
    close($fh);
}


sub rc($){
    my $seq = (reverse(shift));
    $seq =~ tr/ACTGRYSWKMBVDHactgryswkmbvdh/TGACtYRWSMKVBHDgacyrwsmkvbhd/;
    return $seq;
}
