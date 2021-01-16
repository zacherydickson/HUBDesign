#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;

if(@ARGV < 1){
    die "Usage: ".basename($0). " Tree > ClustInfo\n";
}

my $IN = Bio::TreeIO->new(-file => shift);
my $TreeObj = $IN->next_tree; 

my $root = $TreeObj->get_root_node;

foreach my $node ($root->get_all_Descendents){
    foreach my $desc ($node->get_all_Descendents){
        print $node->id,"\t",$desc->id,"\n" if($desc->is_Leaf);
    }
}
