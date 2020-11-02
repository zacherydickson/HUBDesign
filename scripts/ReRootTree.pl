#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;

if(@ARGV < 1){
    die "Usage: ".basename($0)." TreeFile OutputFile\n";
}

my($inFile,$outFile) = @ARGV;

die "In and Out File cannot be the same\n" if($inFile eq $outFile);

my $IN = Bio::TreeIO->new(-file => shift(@ARGV));
my $OUT = Bio::TreeIO->new(-file => ">".shift(@ARGV), -format => 'newick');


foreach my $treeObj ($IN->next_tree){
    my $root = $treeObj->get_root_node;
    my @children = $root->each_Descendent;
    if(scalar(@children) >= 2){ #Tree is only the root and up to one child, no other way to organize
        @children = sort {$b->branch_length <=> $a->branch_length} @children;
        my $CurrentMaxBranch = $children[0]->branch_length + $children[1]->branch_length;
        my %Max = (value => -1, node => undef);
        foreach my $node ($treeObj->get_nodes){
            my $len = $node->branch_length;
            $len = 0 unless(defined $len);
            @Max{qw(value node)} = ($len, $node) if($len > $Max{value});
        }
        unless($Max{value} <= $CurrentMaxBranch){
            $treeObj->reroot_at_midpoint($Max{node},undef);
            $treeObj->splice($root);
        }
    }
    $OUT->write_tree($treeObj);
}
