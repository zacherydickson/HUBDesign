#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Getopt::Std;
use Class::Struct;

my %opts;
getopts('p:o:uh',\%opts);

if (@ARGV < 2 or exists $opts{h}){
    my $Usage = basename($0)." [-options] ClustAssignment ClustSeq1 ... > POSeq";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of txid assignments for consensus sequences for gene clusters,\n".
            "\t\tConstructs  PSeudoGeneomes for each txid\n".
            "===Input\n".
            "\tClustAssignment\tA Tab delimited file with 3 columns: ClusterID, TxID, and penetrance\n".
            "\t\tWhere penetrance is the % of all descendents of txid which contain the cluster\n".
            "\tClustSeq1 ...\tA set of fasta files which have ids which appear in the ClustID column of ClustAssignment\n".
            "===Output\n".
            "\tA Fasta formatted file with a sequence for each unique txid in ClustAssignment.\n".
            "\t\tIndividual Clusters are concatenated together with a sing 'N' separating them\n".
            "===Options\n".
           "\t-p (0-100)\tA Floating point value indicating the minimum penetrance for a cluster [Default: 50.0]\n".
           "\t-o PATH\tThe file to which to output Pseudo Organism Information [Default: POInfo.tab]\n".
           "\t\tPOInfo is a tab delimited file with 4 columns txid clustid cluststart clustlen\n".
           "===Flags\n".
           "-u\tCreate a pseudo-organsism from all unassigned sequences\n".
           "-h\tDisplay this message and exit\n";
    } else {
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

$opts{p} = 50 unless(exists $opts{p});
$opts{o} = "POInfo.tab" unless(exists $opts{o});
open (my $ofh, "| sort > $opts{o}") or die "Could not write to $opts{o}: $!\n";

my ($ClustAssignFile,@SeqFileList) = @ARGV;

##Load ClusterAssignments
my %ClustAssign; #Hash with keys of ClustID and values of txid
if(open(my $fh, $ClustAssignFile)){
    while(my $line = <$fh>){
        chomp($line);
        my ($clustID,$txid,$penetrance) = split(/\t/,$line);
        if($penetrance > $opts{p}){
            $ClustAssign{$clustID} = $txid;
        } elsif(exists $opts{u}){
            $ClustAssign{$clustID} = "unassigned";
        }
    }
    close($fh);
} else {
    die "Could not open $ClustAssignFile: $!\n";
}

##Load and Process Sequences
my %PseudoOrgs; #Hash with keys of txid and values of hash ref with keys of pos and seq, pos has a value of integer and seq has a value of array ref with elements of strings;
foreach my $file (@SeqFileList){
    my $IN = Bio::SeqIO->new(-file => $file, -format => 'fasta');
    while(my $seqObj = $IN->next_seq){
        next unless(exists $ClustAssign{$seqObj->id});
        my $txid = $ClustAssign{$seqObj->id};
        $PseudoOrgs{$txid} = {pos => 1,seq => []} unless(exists $PseudoOrgs{$txid});
        push(@{$PseudoOrgs{$txid}->{seq}},$seqObj->seq);
        print $ofh join("\t",($txid,$seqObj->id,$seqObj->length,$PseudoOrgs{$txid}->{pos})),"\n";
        $PseudoOrgs{$txid}->{pos} += $seqObj->length + 1;
    }
}
close($ofh);

##Output PseudoOrgs
my $OUT = Bio::SeqIO->newFh(-format => 'fasta');
foreach my $txid (sort keys %PseudoOrgs){
    my $seqObj = Bio::Seq->new(-id => $txid, -seq => join("N",@{$PseudoOrgs{$txid}->{seq}}));
    print $OUT $seqObj;
}

