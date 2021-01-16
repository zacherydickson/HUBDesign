#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;
use Class::Struct;

struct (CLUSTER => {id => '$', length => '$', pos => '$'});

sub DetermineCluster($$);
sub GetTxIDList($);

if(@ARGV < 3){
    die "Usage: ".basename($0). " POInfo RegionData Tree > InfoFile\n";
}

my ($POInfoFile, $RegionInfoFile, $TreeFile) = @ARGV;

##LoadPOInfo (4 tab delim columns: txid, clustid, clustlen, and pos within PO
my %POInfo; # Hash with keys of txid and values of array ref with elements of clusterObj
if(open(my $fh,$POInfoFile)){
    while(my $line = <$fh>){
        chomp($line);
        my($tID,$cID,$len,$pos) = split(/\t/,$line);
        $POInfo{$tID} = [] unless(exists $POInfo{$tID});
        push(@{$POInfo{$tID}},CLUSTER->new(id => $cID, length => $len, pos => $pos));
        @{$POInfo{$tID}} = sort {$a->pos <=> $b->pos} @{$POInfo{$tID}};
    }
    close($fh);
} else {
    die "Could not open $POInfoFile: $!\n";
}



##Load Tree
my $IN = Bio::TreeIO->new(-file => $TreeFile, -format => 'newick');
my $TreeObj = $IN->next_tree;

#Desired Output 5 Tab delim columns: regID, txidList, clustID, pos in cluster, lcaTxID, pos in PO

##Process RegionInfo (4 tab delim columns: regID, txID, pos within PO, seq (not used)
if(open(my $fh,$RegionInfoFile)){
    while(my $line = <$fh>){
        chomp($line);
        my($rID,$tID,$posinPO) = split(/\t/,$line);
        my ($cID,$posinClust) = DetermineCluster($posinPO,$POInfo{$tID});
        unless(defined $cID and defined $posinClust){
            warn "[WARNING] Could not identify cluster of $rID: Skipping\n";
            next;
        }
        my $txidStr = GetTxIDList($tID);
        print join("\t",($rID,$txidStr,$cID,$posinClust,$tID,$posinPO)),"\n";
    }
    close($fh);
} else {
    die "Could not open $RegionInfoFile: $!\n";
}


sub DetermineCluster($$){
    my ($pos, $clustListRef) = @_;
    for(my $i = 0; $i < @$clustListRef; $i++){
        my $clustObj = $clustListRef->[$i];
        if($pos >= $clustObj->pos and $pos <= $clustObj->pos + $clustObj->length - 1){
            return ($clustObj->id,$pos - $clustObj->pos + 1);
        }
    }
    return (undef,undef);
}


sub GetTxIDList($){
    my ($txid) = @_;
    my ($node) = $TreeObj->find_node(-id => $txid);
    print STDERR "$txid\n" unless(defined $node);
    if($node->is_Leaf){
        return $txid;
    } else {
        my @txidList;
        foreach my $desc ($node->get_all_Descendents){
            if($desc->is_Leaf){
                push(@txidList,$desc->id);
            }
        }
        return join(",",@txidList);
    }
}
