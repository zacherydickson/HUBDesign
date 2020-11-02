#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;
use Getopt::Std;

sub PartitionTipLabels($$);

my %opts;
getopts('d:',\%opts);

$opts{d} = (exists $opts{d}) ? $opts{d} : 15;

if(@ARGV < 1){
    die "Usage: @{[basename($0)]} -d [0-100] tree\n\t-d\tMax root-tip distance within a cluster [Default: 15%]\n";
}


my $IN = Bio::TreeIO->new(-file => shift(@ARGV), -format => 'newick');

my $treeObj = $IN->next_tree;


my @ClusteredLabels = PartitionTipLabels($treeObj->get_root_node,$opts{d});

my $ClustID = 0;
my %LabelDict;
foreach my $cluster (@ClusteredLabels){
    foreach my $Label (sort {$a <=> $b} @{$cluster}){
        $LabelDict{$Label} = $ClustID;
    }
    $ClustID++;
}

print "$_\t$LabelDict{$_}\n" foreach(sort {$a <=> $b} keys %LabelDict);

sub PartitionTipLabels($$){
    my ($node,$thresh) = @_;
    if($node->is_Leaf){
        return [$node->id];
    }
    if($node->height <= $thresh){
        my @leaves = $node->get_all_Descendents;
        for(my $i = 0; $i < @leaves; $i++){
            unless($leaves[$i]->is_Leaf){
                splice(@leaves,$i,1);
                $i--;
            }
        }
        return ([map {$_->id} @leaves]);
    }
    return map {PartitionTipLabels($_,$thresh)} $node->each_Descendent;
}
