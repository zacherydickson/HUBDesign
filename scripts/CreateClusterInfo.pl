#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

if(@ARGV < 2){
    die "Usage: ".basename($0)." GeneInfo Clust1 ... > ClusterInfo.tab\n".
    "\tNote: It is assumed that seqnumber is the same as in the order of sequences for each gene in GeneInfo, and clust files aare named gene_name.clust, Clsuter Ids are assumed to be numeric\n";
}

my $GeneInfoFile = shift(@ARGV);
my %TxIDinGene; #Hash with keys of Gene Name and values of array ref with elements of txid
if(open(my $fh, $GeneInfoFile)){
    while(my $line = <$fh>){
        chomp($line);
        my($name,$txid) = split(/\t/,$line); 
        $TxIDinGene{$name} = [] unless(exists $TxIDinGene{$name});
        push(@{$TxIDinGene{$name}},$txid);
    }
    close($fh);
} else {
    die "Could not open $GeneInfoFile: $!\n";
}

foreach my $file(@ARGV){
    my $name = basename($file,".clust");
    next unless exists($TxIDinGene{$name});
    if(open(my $fh, $file)){
        my %TxidinCluster;
        while(my $line = <$fh>){
            chomp($line);
            my($seqIndex,$ClustID) = split(/\t/,$line);
            $TxidinCluster{$ClustID} = [] unless(exists $TxidinCluster{$ClustID});
            push(@{$TxidinCluster{$ClustID}},$TxIDinGene{$name}->[$seqIndex-1]);
        }
        close($fh);
        foreach my $ClustNo (sort {$a <=> $b} keys %TxidinCluster){
            foreach my $txid (@{$TxidinCluster{$ClustNo}}){
                print "${name}_".sprintf("%04d",$ClustNo),"\t$txid\n";
            }
        }
    } else {
        warn "Could not open $file: $!; Skipping\n";
    }
}
