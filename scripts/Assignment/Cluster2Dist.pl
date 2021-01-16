#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;

sub initProgress;
sub updateProgress;

my %opts;
getopts('p:v',\%opts);

my @Priority = (1) x @ARGV;
@Priority = split(/,/,$opts{p}) if(exists $opts{p});
my $bVerbose = (exists $opts{v}) ? 1 : 0;


die "Number of Priorities must be the same as the number of files provided\n" if(@Priority < @ARGV);

if(@ARGV < 1){
    die "Usage: ".basename($0). " [-p Priority1,... -v(erbose)] ClusterInfo1 ... > dist.phylip\n".
    "\tPriority Is the Weight[0-1] each set of cluster information is given (Default equal weight)\n".
    "\tAssumption: All Cluster IDs are unique across all files\n";
}

## Read ClusterInfo

my %Cluster;
my %ClusterPriority;
my %txidCount;

while(my $clustFile = shift(@ARGV)){
    print STDERR "Reading $clustFile...\n" if($bVerbose);
    my $priority = shift(@Priority);
    if(open(my $fh,$clustFile)){
        while(my $line = <$fh>){
            chomp($line);
            my($clustID,$txid) = split(/\t/,$line);
            $Cluster{$clustID} = [] unless(exists $Cluster{$clustID});
            if(exists $ClusterPriority{$clustID} and $ClusterPriority{$clustID} != $priority){
                warn "$clustID Observed in multiple files: Assigning priority=$priority\n";
            }
            $ClusterPriority{$clustID} = $priority;
            $txidCount{$txid} = 0 unless(exists $txidCount{$txid});
            $txidCount{$txid}++;
            push(@{$Cluster{$clustID}},$txid);
        }
        close($fh);
    } else {
        die "Could not open $clustFile: $!\n";
    }
}


print STDERR "Allocating Distance Matrix...\n" if($bVerbose);

##Allocate DistMat
my @DistMat;
my @txidList = sort keys %txidCount;
my %index;
@index{@txidList} = (0 .. $#txidList);
foreach my $i (1 .. scalar(@txidList)){
    my @row = (0) x scalar(@txidList);
    push(@DistMat,\@row);
}


print STDERR "Counting Co-occurences...\n" if($bVerbose);

##Get Co-occurence Counts
my $clustCount = scalar(keys %Cluster);
my $curCluster = 0;
my $UpdateAt = $clustCount / 100;
initProgress() if($bVerbose);
while( my ($clustid,$cList) = each %Cluster){
    if($bVerbose and $curCluster++ > $UpdateAt){
        updateProgress;
        $UpdateAt += $clustCount / 100;
    }
    for(my $i = 0; $i < @{$cList}; $i++){
        my $iTxID = $cList->[$i];
        for(my $j = $i + 1; $j < @{$cList}; $j++){
            my $jTxID = $cList->[$j];
            $DistMat[$index{$iTxID}]->[$index{$jTxID}] += $ClusterPriority{$clustid};
        }
    }
}


print STDERR "\nCalculating Distance...\n" if($bVerbose);

##Convert Co-occurence count to Distance
#(Note: Distance = 1 - Co-occurence / min (seqcount))
for (my $i = 0; $i < @txidList; $i++){
    for (my $j = 0; $j < @txidList; $j++){
        if($i == $j){
            $DistMat[$i]->[$j] = 0;
            next;
        }
        if($j < $i){
            $DistMat[$i]->[$j] = $DistMat[$j]->[$i];
            next;
        }
        my $minCount = $txidCount{$txidList[$i]};
        $minCount = $txidCount{$txidList[$j]} if($txidCount{$txidList[$j]} < $minCount);
        warn "[WARN] ZeroCount at $i,$j\n" if($minCount == 0);
        $DistMat[$i]->[$j] /= $minCount if($minCount != 0); 
        $DistMat[$i]->[$j] = 1 - $DistMat[$i]->[$j];
        $DistMat[$i]->[$j] = 0 if($DistMat[$i]->[$j] < 0);
    }
}


print STDERR "Outputing Distances\n" if($bVerbose);

##Output Dist Mat 
print scalar(@txidList),"\n";
for (my $i = 0; $i < @txidList; $i++){
    printf("%-12s%s\n",$txidList[$i],join(" ",@{$DistMat[$i]}));
}


sub initProgress{
    print STDERR join("",((' ') x 101,"]\r[>"));
}
sub updateProgress{
    print STDERR "\b=>";
}
