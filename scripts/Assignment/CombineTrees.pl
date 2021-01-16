#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;

sub removeSinglets($);
sub GetLeafStr($);
sub GetLeafStrDict($);

if(@ARGV < 2){
    die "Usage: ".basename($0). " ClusterTree LineageTree OutputTree\n";
}

my ($clust,$lineage,$out) = @ARGV;

my $IN = Bio::TreeIO->new(-file => $clust);
my $clustTree = $IN->next_tree;
$IN = Bio::TreeIO->new(-file => $lineage);
my $lineageTree = $IN->next_tree;

my $ZeroPad = int(log($lineageTree->number_nodes) / log(10)) + 1;

my $nextID = 0;
foreach my $node ($clustTree->get_root_node->get_all_Descendents){
    if($node->branch_length == 0.5){
        $node->ancestor($clustTree->get_root_node);
    }
    $node->branch_length(undef);
    $node->id(sprintf("TMP%0*d",$ZeroPad,$nextID++)) unless(defined $node->id);
}

removeSinglets($clustTree->get_root_node);

my %lineageLeafStrDict = GetLeafStrDict($lineageTree);
#print GetLeafStr($clustTree->get_root_node),"\n";
#print GetLeafStr($lineageTree->get_root_node),"\n";
foreach my $node ($clustTree->get_nodes){
    next if($node->is_Leaf);
    my $leafStr = GetLeafStr($node);
    my $lineageLabel = $lineageLeafStrDict{$leafStr};
    $node->id($lineageLabel) if(defined $lineageLabel);
}

my $OUT = Bio::TreeIO->new(-file => ">$out");
$OUT->write_tree($clustTree);

sub removeSinglets($){
    my $node = shift;
    return if($node->is_Leaf);
    my @children = $node->each_Descendent;
    if(scalar(@children) == 1){
        my ($child) = @children;
        $child->ancestor($node->ancestor);
        $node->ancestor->remove_Descendent($node);
        removeSinglets($child);
    } else {
        removeSinglets($_) foreach (@children);
    }
    
}

sub GetLeafStr($){
    my $node = shift;
    my @leafLabelList;
    foreach my $desc ($node->get_all_Descendents){
        push(@leafLabelList,$desc->id) if($desc->is_Leaf);
    }
    return join(";",(sort @leafLabelList));
}

sub GetLeafStrDict($){
    my $treeObj = shift;
    my %LeafStrDict;
    foreach my $node ($treeObj->get_nodes){
        next if($node->is_Leaf);
        $LeafStrDict{GetLeafStr($node)} = $node->id;
    }
    return %LeafStrDict;
}
