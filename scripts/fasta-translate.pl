#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

sub ProcessNumericOption($$$$$$);

my %opts;
getopts('c:f:t:h',\%opts);

if(@ARGV < 1 or exists $opts{h}){
    my $Usage = basename($0). " [-options] DNA1.fna[.gz] ... > Prot.fna";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of DNA sequences describing mRNA, outputs the translated protein sequences\n".
            "===Usage\n\t$Usage\n".
            "===Inputs\n".
            "\tDNA#.fna[.gz]\tA fasta formatted file (optionally gzipped) containing DNA sequences describing mRNA\n".
            "===Output\n".
            "\tA Single Fasta Formated file containing the translations from all input files\n".
            "===Options\n".
            "\t-c INT\tAn Integer Describing which genetic code to use (Default: 1)\n".
            "\t-f INT\tAn Integer describer the frame to translate in (Default: 0)\n".
            "\t-t CHAR\tThe Character to be used as a stop codon (Default: *)\n".
            "===Flags\n".
            "\t-h\tDisplay this message and exit\n";
    } else {
        die "Usage: $Usage\nUse -h for more info\n";
    }
}

$opts{c} = ProcessNumericOption($opts{c},1,1,33,1,"Genetic Code");
$opts{f} = ProcessNumericOption($opts{f},0,undef,undef,1,"Frame");
$opts{t} = (exists $opts{t}) ? substr($opts{t},0,1) : '*';

my @FileList = @ARGV;

my $Out = Bio::SeqIO->newFh(-format => 'fasta');
foreach my $file (@FileList){
    my $opencmd = ($file =~ m/\.gz$/) ? "gzip -dc $file |" : $file;
    open (my $fh, $opencmd) or die "Could not open $file: $!\n";
    my $IN = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    while(my $seqObj = $IN->next_seq){
        print $Out $seqObj->translate(-codontable_id => $opts{c}, -frame => $opts{f}, -terminator => $opts{t})
    }
    close($fh);
}

sub ProcessNumericOption($$$$$$){
    my($val,$default,$min,$max,$bInt,$varName) = @_;
    return $default unless(defined $val);
    if(looks_like_number($val)){
        if(!$bInt or $val == int($val)){
            if(!defined $min or $val >= $min){
                if(!defined $max or $val <= $max){
                    return $val;
                }
            }
        }
    }
    my $message = sprintf("%s must be a%s between %s and %s",$varName,
       ($bInt ? "n integer" : " value"),(defined $min ? $min : "-∞"),(defined $max ? $max : "∞")); 
    die "[ERROR] $message\n";
}
