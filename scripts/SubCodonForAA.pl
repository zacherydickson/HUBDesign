#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;

if (@ARGV < 2) {
    die "Usage: ".basename($0)." CodonSeq AlignedProtSeq\n";  
}

my($CodonFile, $AlnFile) = @ARGV;

my $AlnIn = new Bio::SeqIO(-file => $AlnFile, -format => 'fasta');

##Load AlignmentInfo
my %AlnDesc; #Hash with keys of seqID and values of array ref, with an odd number elements alternating between a gap length and a sequence length; ends and begins with a gap
while(my $seqObj = $AlnIn->next_seq){
    my $id = $seqObj->id;
    $id =~ s/_1$//;
    my $seq = $seqObj->seq;
    my $bGap = 1;
    my $run = "";
    while(my $char = substr($seq,0,1,"")){
        my $res = ($char eq "-");
        $res = !$res unless($bGap);
        if($res){
            $run .= $char;
        } else{
            push(@{$AlnDesc{$id}},length($run));
            $run = $char;
            $bGap = !$bGap;
        }
    }
    push(@{$AlnDesc{$id}},length($run));
    push(@{$AlnDesc{$id}},0) unless($bGap);
}

#while(my ($id, $arrRef) = each %AlnDesc){
#    print ">$id\n";
#    print join(" ",@{$arrRef}),"\n";
#}
#exit;
##Process Codons
my $CodonIn = new Bio::SeqIO(-file => $CodonFile, -format => 'fasta');
my $CodonOut = Bio::SeqIO->newFh(-format => 'fasta');
while( my $seqObj = $CodonIn->next_seq){
    next unless(exists $AlnDesc{$seqObj->id});
    my $bGap = 1;
    my $seq = $seqObj->seq;
    my $codonaln = "";
    while(@{$AlnDesc{$seqObj->id}}){
        my $runLen = shift(@{$AlnDesc{$seqObj->id}});
        if($bGap){
            $codonaln .= "-" x ($runLen * 3);
        } else {
            my $cdnLen = $runLen * 3;
            $cdnLen = length($seq) if $cdnLen > length($seq);
            $codonaln .= substr($seq,0,$runLen * 3,"");
        }
        $bGap = !$bGap;
    }
    my $stopLen = (length($seq) >= 3) ? 3 : length($seq);
    $codonaln .= substr($seq,0,$stopLen,"");
    if(length($seq)){
        warn "[WARNING] All of sequence @{[$seqObj->id]} not used: Remaining $seq\n";
    }
    my $outSeqObj = new Bio::Seq(-seq => $codonaln, -id => $seqObj->id);
    print $CodonOut $outSeqObj;
}
