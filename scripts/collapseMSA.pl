#!usr/perl
use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;

sub val_in_array($@);

my $Usage = "Usage: perl collapseMSA.pl [-options] msa.fasta >out\n".
            "\t -l [Int] : Minimum unique segment length\n".
            "\t -c : Output simple consensus sequence\n".
            "\t -o : Output conserved regions\n".
            "\t -a [String] : The Header to display for collapsed sequence\n".
            "\t -g [d|k|fd|fk]: Consensus Gap behaviour: d=drop, k=keep, f=fall_through\n".
            "\t\t fall_through means to go to the next most common character [Default: d]\n".
            "\t -f [float]: Use with '-g f'; When falling through, the minimum proportion of sequences\n".
            "\t\twhich agree on a non-gap, otherwise the gap is dropped/kept (0 - 0.5) [Default: 0]\n".
            "\t -m [n,g,p]: The Consensus Value to report n=none, g=geo_mean of support for majority nucleotide\n".
            "\t\t p=proportion of divergent positions [Default n]\n".
            "\t -h : Output this message\n";


my %opts;
getopts('l:ca:o:g:hf:m:',\%opts);

die $Usage if (@ARGV < 1 or exists $opts{h});

my $msaFile = shift(@ARGV);

my $oLength = (exists $opts{l}) ? $opts{l} : 75;
my $bConsensus = (exists $opts{c}) ? 1 : 0;
my $SharedHeader =  (exists $opts{a}) ? $opts{a}: "Collapsed";
my $bGapDrop = 1;
my $bGapKeep = 0;
my $bSecondaryDrop = 1;
my $bReportCons = 0;
my $bReportGeoMean = 1;
if (exists $opts{g}){
    my %tmp;
    @tmp{qw(d k fk fd)} = (1) x 4;
    unless(exists $tmp{$opts{g}}){
        die "Invalid -g option $opts{g}; must be one of d,k,fd, or fk\n";
    }
    if($opts{g} ne 'd'){
        $bGapKeep = 1 if($opts{g} eq 'k');
        $bGapDrop = 0;
        $bSecondaryDrop = 0 if($opts{g} eq 'fk');
    }
}
if(exists $opts{m}){
    my %valid;
    @valid{qw(n g p)} = (1) x 3;
    unless(exists $valid{$opts{m}}){
        die "Invalid -m option $opts{m}; must be one of c g or n\n";
    }
    $bReportCons = 1 if($opts{m} ne 'n');
    $bReportGeoMean = 0 if($opts{m} eq 'p');
}
my $MinFallProp = (exists $opts{f}) ? $opts{f} : 0;
if($MinFallProp > 0.5){
    print STDERR "[WARNING] -f must be less than or equal to 0.5, continuing with 0.5\n";
    $MinFallProp = 0.5;
}


unless($bConsensus){
    die "Non-naive consensus is not yet implemented, please use -c or use a different program\n";
}

my @headerList;
my $seq;
my @charMatrix;

my $msaLength;
my $in = Bio::SeqIO->new(-file => $msaFile, -format => 'fasta');
while(my $seqObj = $in->next_seq){
    my $header = $seqObj->id."|".$seqObj->desc;
    push(@headerList,$header);
    my @tmp = split('',$seqObj->seq);
    if(defined($msaLength)){
        unless(scalar(@tmp) == $msaLength){
            die "Length of " . $headerList[$#headerList] . " does not match previous lengths in $msaFile\n";
        }
    } else{
        $msaLength = scalar(@tmp);
    }
    push(@charMatrix,\@tmp);
}

unless(defined $msaLength){
    die "[ERROR] No sequences found in $msaFile\n";
}

my $consensus;
my $conservation = 0;
my $divergeProp = 0;
my $seqCount = scalar(@headerList);
my @commonList = ('_') x $msaLength;
my $bFirst = 0;
for (my $j = 0; $j < $msaLength; $j++){ 
    my %charCount;
    my @seqList;
    my $max = -1;
    my $maxchar;
    my $undermaxchar;
    my $undermax = -1;
    for (my $i = 0; $i < $seqCount; $i++){
        my $char = $charMatrix[$i][$j];
        if(exists $charCount{$char}){
            $charCount{$char}++;
        } else {
            $charCount{$char}=1;
            push(@seqList,$i);
        }
    }
    foreach my $char (sort keys %charCount){
        if($charCount{$char} >= $max){
            $max = $charCount{$char};
            $maxchar = $char;
        }
    }
    foreach my $char (sort keys %charCount){
        if($char ne $maxchar and $charCount{$char} >= $undermax){
            $undermax = $charCount{$char};
            $undermaxchar = $char;
        }
    }
    #my @arr = keys %charCount;
    #print STDERR "@arr\n" if(@arr > 1);
    #print STDERR "$j\n" if(join('',@arr) eq '-A');
    #print STDERR "$j: $maxchar) $max, $undermax, $seqCount, $maxchar\n";
    my $nonGapProp;
    $nonGapProp = $undermax/($seqCount) if($maxchar eq '-');
    my $gapCount = (exists $charCount{'-'} ? $charCount{'-'} : 0);
    if($maxchar ne '-' or $bGapKeep or ($nonGapProp < $MinFallProp and !$bSecondaryDrop)){
    	$consensus .= $maxchar;
        $divergeProp++ if($max < $seqCount - $gapCount);
	$conservation += log($max / $seqCount);
    } elsif(!$bGapDrop and $nonGapProp > $MinFallProp) {
        $consensus .= $undermaxchar;
        $divergeProp++ if($max < $seqCount - $gapCount);
	$conservation += log($undermax/$seqCount);
    }# elsif ($max < 2 and $seqCount > 1){
    #    print STDERR "$SharedHeader\t$j\n";
    #}
    my @tmp = keys(%charCount);
    if(@tmp == 1){
        $commonList[$j] = $tmp[0];
    } else {
        $commonList[$j] = \@seqList;
    }
    #if(!$bFirst){
    #    print STDERR "$j:$conservation\n";
    #}
    #if($conservation == 0 and !$bFirst++){
    #    print STDERR "$j: $maxchar-$max, $undermaxchar-$undermax\n";
    #}
}
$conservation = exp($conservation / length($consensus)) * 100;
$divergeProp /= length($consensus);
$divergeProp *= 100;

if($bConsensus){
    my $out = Bio::SeqIO->newFh(-format => 'fasta');
    my $header = $SharedHeader;
    if($bReportCons){
        $header .= ";".($bReportGeoMean ? $conservation : $divergeProp);
    }
    my $seqObj = Bio::Seq->new(-id => $header, -seq => $consensus);
    print $out $seqObj;
    #print ">$SharedHeader-$conservation\n";
    #print "$consensus\n";
} else {
    die "Non Consensus use is not yet implemented, please use -c\n";
    my ($seq2, $tmpseq);
    my $lastLen = 0;
    my $laststart = 0;
    my $state = 0;
    for(my $j= 0; $j < scalar(@commonList); $j++){
        my $tmpstate = (ref($commonList[$j]) eq "ARRAY" and val_in_array(2,@{$commonList[$j]}));
        if($tmpstate != $state and defined($tmpseq)){
            if($state == 1){ #From unique to common
                if(length($tmpseq) >= $oLength){
                    $lastLen = length($tmpseq);
                    $seq2 .= $tmpseq;
                    $tmpseq = '';
                    $state = 0;
                }
            } else { #From common to unique
                if(length($tmpseq) >= $oLength){
                    $tmpseq = '';
                    $seq2 .= "\n>seq2\n";
                    $laststart = $j;
                } else {
                    $tmpseq = substr($seq2,$laststart,$lastLen,'').$tmpseq;
                }
                $state = 1;
            }
        } elsif (! defined($tmpseq)){
            $state = $tmpstate;
        }
        $tmpseq .= $charMatrix[2][$j] if($charMatrix[2][$j] ne '-');
    }
    if($state == 1){
        $seq2 .= $tmpseq;
    }
    
    print ">$SharedHeader\n";
    for(my $j= 0; $j < scalar(@commonList); $j++){
       if(ref($commonList[$j])){
           if(val_in_array(2,@{$commonList[$j]})) {
               print "_";
           } else {
               print $charMatrix[0][$j];
           }
       } else{
           print $commonList[$j];
       }
    }
    print "\n>seq2\n$seq2";
    print "\n";
}

sub val_in_array($@){
    my ($val,@arr) = @_;
    my %hash;
    @hash{@arr} = (0) x scalar(@arr);
    return exists($hash{$val});
}
