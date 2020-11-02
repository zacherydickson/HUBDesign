#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::Seq;

if (@ARGV < 1){
    die "Usage: ".basename($0)." reshapedOligFile\n";
}

sub rc{
    return Bio::Seq->new(-seq => shift, -id => 'arb')->revcom()->seq;
}

my %Struct;
my %TgtDict;
if(open(my $fh, shift(@ARGV))){
    while(my $line = <$fh>){
        chomp($line);
        my($id,$tgt,$pos,$seq) = split(/\t/,$line);
        my $rcseq = rc($seq);
        unless(exists $Struct{$seq} || exists $Struct{$rcseq}){
            $Struct{$seq} = [];
        }
        if(exists $Struct{$rcseq}){
            $seq = $rcseq;
        } 
        $TgtDict{$id} = $tgt;
        push(@{$Struct{$seq}},$id);
    }
} else {
    die "Could not open input: $!\n";
}


#Each Row of output is baitID tab DupID
my $DupID = 0;
foreach my $seq (sort {scalar(@{$Struct{$b}}) <=> scalar(@{$Struct{$a}})} keys %Struct){
    last if(scalar(@{$Struct{$seq}}) <= 1);
    $DupID++;
    my %UniqTgts;
    @UniqTgts{(map {$TgtDict{$_}} @{$Struct{$seq}})} = (1) x @{$Struct{$seq}};
    if(scalar(keys %UniqTgts) > 1){
        warn "[WARNING] Duplicates found across targets: $DupID hits @{[keys %UniqTgts]}\n";
    }
    foreach my $baitID (@{$Struct{$seq}}){
    print $baitID,"\t",sprintf("Dup%04d",$DupID),"\n";
    }
}

