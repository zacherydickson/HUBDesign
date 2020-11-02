#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::TreeIO;
use Getopt::Std;
use Parallel::ForkManager;
use Scalar::Util qw(looks_like_number);

sub unique(@);
sub get_leaf_descendents($);
sub GetProcessorCount();
sub ProcessNumericOption($$$$$$);

my %opts;
getopts('t:vh',\%opts);

my $max_proc = ProcessNumericOption($opts{t},1,0,undef,1,"Threads");
my $cpu_count = GetProcessorCount();
if($cpu_count == -1){
    warn "[WARNING] Could not determine sys_max threads: $!\n".
    "\t Can't fully validate max_threads\n";
    if($max_proc == 0){
        warn "[WARNING] sys_max threads unknown: proceeding with 1 thread\n";
        $max_proc = 1;
    }
} elsif($max_proc == 0 or $max_proc > $cpu_count){
    if($max_proc > $cpu_count){
        warn "[WARNING] Max threads is greater than sys_max: proceeding with $cpu_count threads\n";
    }
    $max_proc = $cpu_count;
}
my $bVerbose = exists $opts{v} ? 1 : 0;

if(@ARGV < 2 or exists $opts{h}){
    die "Usage: ".basename($0)."[-options] GeneInfo.tab Lineage.tree\n".
        "\t-t INT\tNumber of Threads; 0 = sys_max [Default 1]\n".
        "\t-v\tVerbose output\n".
        "\t-h\tDisplay help and exit\n";
}

my($GeneFile,$LineageFile) = @ARGV;

##Load Gene Info
print STDERR "[INFO] Loading Gene Info\n" if($bVerbose);
my %TxIDSetDict; # Hash with keys of gene and values of array ref with elements of txid's
if(open(my $fh, $GeneFile)){
    while(my $line = <$fh>){
        chomp($line);
        next if($line =~ /^$/);
        my($Gene, $txid) = split(/\t/,$line);
        $TxIDSetDict{$Gene} = [] unless(exists $TxIDSetDict{$Gene});
        push(@{$TxIDSetDict{$Gene}},$txid);
    }
    close($fh);
} else {
    die "Could not open $GeneFile: $!\n";
}

##Load Lineage Tree
print STDERR "[INFO] Loading Lineage\n" if($bVerbose);
my $IN = Bio::TreeIO->new(-file => $LineageFile, -format => 'newick');
my $treeObj = $IN->next_tree;


my $pm = Parallel::ForkManager->new($max_proc);
    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $cor_dump, $data_structure_reference) = @_;
            if($exit_code == 0){
                if (defined $data_structure_reference){
                    print "$_\n" foreach(@{$data_structure_reference});
                    print STDERR "[INFO] The block beginning at $ident has been completed\n" if($bVerbose);
                } else {
                    warn "[WARNING] No result from $ident\n";
                }
            } else {
                warn "[WARNING] Unexpected result from block beginning at $ident\n";
            }
        }
    );

##Process Genes
my $counter = 0;
my %lcaDict;
my @geneList = keys(%TxIDSetDict);
my $geneCount = scalar(keys %TxIDSetDict);
my $blockSize = int($geneCount / $max_proc);
$blockSize = 1 if ($blockSize < 1);
my $targetUpdateTime = 10 * $max_proc;
GENE:
for (my $blockStart = 0; $blockStart < @geneList; $blockStart += $blockSize){
    $pm->start($blockStart) and next GENE;
    my $blockEnd = $blockStart + $blockSize-1;
    $blockEnd = ($blockEnd > $#geneList) ? $#geneList : $blockEnd;
    my @genesInBlock = @geneList[$blockStart .. $blockEnd];
    my @Output;
    my $lastUpdate = 0;
    my $count = 0;
    my $updateRate = 10; #Per thousand progress per update
    #To introduce an offset ofr each block to make reporting come consistently rather than in blocks
    my $initUpdateOffset = $targetUpdateTime - $updateRate * int($blockStart/$blockSize);
    my $lasttime = time;
    foreach my $gene (@genesInBlock){
        $count++;
        my @txidList = unique(@{$TxIDSetDict{$gene}});
        my $txidKey = join(';',@txidList);
        my $lcaID = $txidList[0];
        my $penetrance = 100;
        if(exists $lcaDict{$txidKey}){
            ($lcaID,$penetrance) = @{$lcaDict{$txidKey}};
        } else {
            my @nodeList;
            push(@nodeList,$treeObj->find_node(-id => $_)) foreach (@txidList); 
            my $leafCount = 1;
            if(@nodeList >= 2){
                my $lcaNode = $treeObj->get_lca(-nodes => \@nodeList);
                $lcaID = $lcaNode->id;
                my @leafList = get_leaf_descendents($lcaNode);
                $leafCount = scalar(@leafList);
                $penetrance = @nodeList/$leafCount*100
            }
            $lcaDict{$txidKey} = [$lcaID,$penetrance];
        }
        my $progress = int(($count / $blockSize) * 1000);
        if($bVerbose and $progress - $lastUpdate >= $updateRate){
            my $curtime = time;
            my $elapsed = $curtime - $lasttime;
            $elapsed = 1 if($elapsed < 1);
            if($elapsed != $targetUpdateTime){
                $updateRate = ($targetUpdateTime - $initUpdateOffset) * ($progress - $lastUpdate) / $elapsed;
                $updateRate = 1 if ($updateRate) < 1;
            }
            $lastUpdate = $progress;
            print STDERR "[INFO] The block beginning at $blockStart is @{[$lastUpdate/10]}% complete (@{[$curtime - $lasttime]}s)\n";
            $lasttime = $curtime;
            $initUpdateOffset = 0;
        }
        my $str = join("\t",($gene,$lcaID,sprintf("%0.2f",$penetrance)));
        push(@Output,$str);
    }
    $pm->finish(0,\@Output);
}
$pm->wait_all_children;

sub unique(@){
    my %Hash;
    @Hash{@_} = (1) x @_;
    return keys %Hash;
}

sub get_leaf_descendents($){
    my $self = shift;
    my @descendents = $self->get_all_Descendents;
    my $count = scalar(@descendents);
    for(my $i = 0; $i < $count; $i++){
        my $node = shift(@descendents);
        push(@descendents,$node) if($node->is_Leaf);
    }
    return @descendents;
}

sub GetProcessorCount(){
    my $cpu_count = -1;
    if(open my $handle, "/proc/cpuinfo"){
        $cpu_count = scalar(map /^processor/, <$handle>);
        close($handle);
    }
    return $cpu_count;
}

sub ProcessNumericOption($$$$$$){
    my($val,$default,$min,$max,$bInt,$varName) = @_;
    return $default unless(defined $val);
    if(looks_like_number($val)){
        if(!$bInt or $val == int($val)){
            if(!defined $min or $val >= $min){
                if(!defined $max or $val <= $max){
                    return $val;
                }
            }
        }
    }
    my $message = sprintf("%s must be a%s between %s and %s",$varName,
       ($bInt ? "n integer" : " value"),(defined $min ? $min : "-∞"),(defined $max ? $max : "∞")); 
    die "[ERROR] $message\n";
}
