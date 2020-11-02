#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;


sub overlap($$$$);

if(@ARGV < 2){
    die "Usage: ".basename($0). " Intervals InInfo\n";
}

my $minRegLen = 75;

my ($IVFile, $InFile) = @ARGV;

##Load Intervals
my %IVDict; #Hash with keys of RegionID and values of array ref with elements of array ref with two elements, start and end
if(open(my $fh, $IVFile)){
    while(my $line = <$fh>){
        chomp($line);
        my ($regid,$s,$e) = split(/\t/,$line);
        $IVDict{$regid} = [] unless(exists $IVDict{$regid});
        ($s,$e) = ($e,$s) unless($s < $e);
        push(@{$IVDict{$regid}},[$s,$e]);
    }
} else {
    die "Could not open $IVFile: $!\n";
}

##Collapse overlapping intervals
while(my ($regid,$ivListRef) = each %IVDict){
    @$ivListRef = sort {$a->[0] <=> $b->[0]} @$ivListRef;
    for(my $i = 1; $i < @$ivListRef; $i++){
        my($s1,$e1) = @{$ivListRef->[$i-1]};
        my($s2,$e2) = @{$ivListRef->[$i-1]};
        my $ol = overlap($s1,$e1,$s2,$e2);
        warn "[WARNING] Interval with end before start in $regid\n" if($ol == -1);
        if($ol){
            $ivListRef->[$i-1] = [$s1,$e2];
            splice(@$ivListRef,$i,1);
            $i--;
        }
    }
}

##Process Regions
if(open(my $fh, $InFile)){
    while(my $line = <$fh>){
        chomp($line);
        my ($regid,$tID,$pos,$seq) = split(/\t/,$line);
        unless(exists $IVDict{$regid}){
            print "$line\n";
            next;
        }
        my @substrArgs; #array of array ref with two elements, start and length
        my $index = 0;
        #Note: ivRef is 1 indexed, substr/$index is 0 indexed
        foreach my $ivRef (@{$IVDict{$regid}}){
            my $len = $ivRef->[0] - $index - 1;
            if($len >= $minRegLen){
                push(@substrArgs,[$index,$len]);
            }
            $index = $ivRef->[1];
        }
        my $remainLen = length($seq) - $index;
        if($remainLen >= $minRegLen){
            push(@substrArgs,[$index,$remainLen]);
        }
        my $subseqcount = 0;
        foreach my $substrarg (@substrArgs){
            $subseqcount++;
            my $newid = "${regid}_$subseqcount";
            my $newpos = $pos+$substrarg->[0];
            my $newseq = substr($seq,$substrarg->[0],$substrarg->[1]);
            print join("\t",($newid,$tID,$newpos,$newseq)),"\n";
        }
    }
} else {
    die "Could not open $IVFile: $!\n";
}

sub overlap($$$$){
    my ($s1,$e1,$s2,$e2) = @_;
    return -1 if($e1 < $s1 or $e2 < $s2);
    my $overlap = 0;
    if($e1 >= $s2 and $s1 <= $e2){
        $overlap = $e1 - $s2 + 1 - (($s1 > $s2) ? ($s1 - $s2) : 0) - (($e1 > $e2) ? ($e1 - $e2) : 0);
    }
    return $overlap;
}
