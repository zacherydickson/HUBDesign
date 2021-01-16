#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV < 2){
    die "Usage: $0 input.fasta output.fasta\n";
}

my $inFile = shift(@ARGV);
my $outFile = shift(@ARGV);

my $In = Bio::SeqIO->new(-file => $inFile, -format => 'fasta');
my $Out = Bio::SeqIO->new(-file => ">$outFile", -format => 'fasta');

my $bSortStrict = 1;
my $baitLen = 75;

my $initPos=-1;
my $curPos=0;
my $curLength=0;
my $curSeq = "";
my $curtxid = "";
my $curtarget = "";
my $adulterated = 0;
while(my $curObj = $In->next_seq){
    my $id = $curObj->id;
    my($txid,$info) = split(/\|/,$id);
    my($targetinfo,%baitData) = split(/#|=/,$info);
    if($txid ne $curtxid or $targetinfo ne $curtarget or $baitData{GenomePos} > $curPos + 1){
        if($initPos != -1){
            my $display_id = "$curtxid|$curtarget#GenomePos=$initPos#Length=$curLength";
            my $seqObj = Bio::Seq->new(-seq => $curSeq, -display_id => $display_id);
            $Out->write_seq($seqObj);
        }
        $initPos = $baitData{GenomePos};
        $curPos = $initPos;
        $curLength = $curObj->length();
        $curSeq = $curObj->seq();
        $curtxid = $txid;
        $curtarget = $targetinfo;
        $adulterated = 0;
    } else {
        if($baitData{GenomePos} < $curPos){
            die "[ERROR] Baits within a txid,target combo are not sorted numerically by Genome Position\n\tEnsure They are sorted properly and try again\n" if($bSortStrict);
            my @TargetNames = split(/:|;/,$targetinfo);
            my %TargetLen = @TargetNames;
            @TargetNames = @TargetNames[map { $_ * 2} 0 .. int($#TargetNames / 2)]; 
            my $TargetPos = $initPos;
            my $i = 0;
            while($TargetPos > $TargetLen{$TargetNames[$i]}){
                $TargetPos -= $TargetLen{$TargetNames[$i]}+1;
                $i++;
            }
            my $SourcePos = $baitData{GenomePos};
            my $j = 0;
            while($SourcePos > $TargetLen{$TargetNames[$j]}){
                $SourcePos -= $TargetLen{$TargetNames[$j]}+1;
                $j++;
            }
            $SourcePos += $baitLen - 1;
            my $insertbase = substr($curObj->seq,-1);
            warn "[WARNING 1] Inserted Base ($insertbase) from $txid|$TargetNames[$i]:$SourcePos at @{[1+$curPos+$adulterated+$baitLen-$initPos]} in $txid|$TargetNames[$i]:$TargetPos\n";
            $adulterated++;
            $curPos--;
        }
        $curPos++;
        $curLength++;
        $curSeq .= substr($curObj->seq,-1);
    }
}
my $display_id = "$curtxid|$curtarget#GenomePos=$initPos#Length=$curLength";
            my $seqObj = Bio::Seq->new(-seq => $curSeq, -display_id => $display_id);
            $Out->write_seq($seqObj);

