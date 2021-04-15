#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use Bio::TreeIO;
use Bio::Matrix::IO;
use Bio::Tree::DistanceFactory;
use IO::String;
use Class::Struct;
use Parallel::ForkManager;

#============================================================
#Declarations

struct  (CLUSTER => {uid => '$', members => '@', assignment => '$', penetrance => '$'});

my $_ID_PADDING = 5;

sub LoadClusterInfo($);      #Usage: my @clusterList = LoadClusterInfo($clustFile);
sub ConstructDendrogram($@); #Usage: my $treeObj = ConstructTaxonTree($guideTreeFile,@clusterList);
sub ReRootTree($);           #Usage: ReRootTree($treeObj);
sub GetProcessorCount();     #Usage: my $sys_proc_count = GetProcessorCount();
sub get_leaf_descendents($); #Usage: my @leafList = get_leaf_descendents($lcaNode);
sub unique(@);               #Usage: my @unique = unique(@list)
sub AssignTaxa($$);          #Usage: AssignTaxa($dendObj,\@clusterList);

#============================================================
#Handling Input

my %opts;
getopts('g:G:w:t:v:h',\%opts);

if(@ARGV < 1 or exists $opts{h}){
    my $Usage =  "Usage: ".basename($0). " [-options] ClusterInfo > Assignment.tsv";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of clusters, the taxa represented by those clusters, constructs\n".
            "\t\ta tree on the taxa, and assigns each cluster to a node in the tree\n".
            "===Usage\n\t$Usage\n".
            "===Input\n".
            "\tClusterInfo\tA tab delimited file with the first two columns being ClusterID and TaxonID\n".
            "===Output\n".
            "\tA tab separated file with columns of ClusterID, TaxonID, and Penetrance\n".
            "===Options\n".
            "-g PATH\tA file containing a newick formatted tree to guide taxon assignment\n".
            "\t\tNote: There is assumed to be no overlap between cluster ids and node ids in the tree\n".
            "-G PATH\tThe file to which to write the dendrogram used to assign clusters to taxa (Default: Guide.dend)\n".
            "-w (0,∞]\tThe relative weight of the topology of the Clustering compared to the\n".
            "\t\tGuide Tree (Default: 1.0)\n".

            "-t INT\tNumber of Threads; 0 = sys_max [Default 1]\n".
            "===Flags\n".
            "-v\tVerbose Output\n";
            "-h\tDisplay this message an exits\n";
    } else{
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

if(exists $opts{g}){
    die "Could not locate file: $opts{g}\n" unless(-e $opts{g});
} else{
    $opts{g} = undef;
}
$opts{G} = (exists $opts{G}) ? $opts{G} : "Guide.dend";
$opts{w} = (exists $opts{w}) ? $opts{w} : 1;
my $max_proc = (exists $opts{t}) ? $opts{t} : 0;
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

#============================================================
#Main Script

my $clustInfoFile = shift(@ARGV);
my @clusterList = LoadClusterInfo($clustInfoFile);
printf STDERR "Clusters: %d\n", scalar(@clusterList) if(exists $opts{v});
print STDERR "\tConstructing dendrogram...\n" if(exists $opts{v});
my $dendObj = ConstructDendrogram($opts{g},@clusterList);
Bio::TreeIO->new(-format => "newick",-file => ">".$opts{G})->write_tree($dendObj);
print STDERR "\tAssigning Taxa...\n" if(exists $opts{v});
AssignTaxa($dendObj,\@clusterList);
print join("\t",qw(ClusterID TaxonID Penetrance)),"\n";
foreach my $clustObj (@clusterList){
    print join("\t",($clustObj->uid,$clustObj->assignment,sprintf("%.2f",$clustObj->penetrance))),"\n";
}
print STDERR "Done\n" if(exists $opts{v});


#============================================================
#Subroutines

#Given a cluster info file, parse out the cluster ids and their associated taxon ids
#Input: a path to a cluster info file
#Output: an array of CLUSTER objects
sub LoadClusterInfo($){
    my ($file) = @_;
    my @clusterList;
    my %indexDict;
    open(my $fh, $file) or die "Could not open ClustInfo file - $file: $!\n";
    my $header = <$fh>;
    while(my $line = <$fh>){
        chomp($line);
        my ($cid, $tid) = split(/\t/,$line);
        if(exists $indexDict{$cid}){
            push(@{$clusterList[$indexDict{$cid}]->members},$tid);
        } else {
            push(@clusterList,CLUSTER->new(uid => $cid, members => [$tid]));
            $indexDict{$cid} = $#clusterList;
        }
    }
    close($fh);
    return(@clusterList);
}

#Given clustering information and a guide tree, constructs a dendrogram 
# with topology in agreement with the observed clustering, with ambiguities resolved by the
# guide tree
#Input: a (possibly undefined) path to a newick formated file, and a list of CLUSTER objects
#Output: a Bio::Tree Object
sub ConstructDendrogram($@){
    my ($guideTreeFile, @clusterList) = @_;
    ##Pass over cluster List to parse out useful info
    my %taxaSet; #The unique taxon ids present only in the cluster info
    my %clusterSet; #The unique cluster ids present only in the cluster info
    foreach my $clustObj (@clusterList){
        $clusterSet{$clustObj->uid} = 1;
        @taxaSet{@{$clustObj->members}} = (1) x scalar(@{$clustObj->members});
    }
    ##Load Information from the Guide Tree
    my %INodeDict; #Hash with keys of ';' delimited taxidLists, and values of txid
    if(defined $guideTreeFile){
        my $treeObj = Bio::TreeIO->new(-file => $guideTreeFile, -format => "newick")->next_tree();
        my $root = $treeObj->get_root_node;
        my $nextID = 0;
        foreach my $node ($root->get_all_Descendents){
            unless(defined $node->id and $node-> id ne ""){
                $node->id(sprintf("TID_%0*d",$_ID_PADDING,$nextID++));
            }
        }
        foreach my $node ($root->get_all_Descendents){
            next if($node->is_Leaf);
            my $clustObj = CLUSTER->new(uid => $node->id);
            foreach my $desc ($node->get_all_Descendents){
                push(@{$clustObj->members},$desc->id) if $desc->is_Leaf;
            }
            $INodeDict{join(';',(sort @{$clustObj->members}))} = $clustObj->uid;
            push(@clusterList,$clustObj);
        }
    } 
    ##Allocate Distance Matrix
    my @DistMat;
    my @txidList = sort keys %taxaSet;
    my %indexDict;
    @indexDict{@txidList} = (0 .. $#txidList);
    foreach my $i (1 .. scalar(@txidList)){
        my @row = (0) x scalar(@txidList);
        push(@DistMat,\@row);
    }
    ##Get Co-occurence Counts, only using taxa fromt the cluster info
    foreach my $clustObj (@clusterList){
        my $weight = (exists $clusterSet{$clustObj->uid}) ? $opts{w} : 1;
        my @tidList = @{$clustObj->members};
        for(my $i = 0; $i < @tidList - 1; $i++){
            my $tid_i = $tidList[$i];
            next unless(exists $taxaSet{$tid_i});
            my $rowIndex = $indexDict{$tidList[$i]};
            for(my $j = $i + 1; $j < @tidList; $j++){
                my $tid_j = $tidList[$j];
                next unless(exists $taxaSet{$tid_j});
                my $colIndex = $indexDict{$tid_j};
                $DistMat[$rowIndex]->[$colIndex] += $weight;
                $DistMat[$colIndex]->[$rowIndex] += $weight;
            }
        }
    }
    ##Determine the number of times each taxon appears
    my %txidCount;
    for(my $i = 0; $i < @txidList; $i++){
        $txidCount{$txidList[$i]} = eval(join("+",@{$DistMat[$i]}));
    }
    ##Convert Co-occurence count to Distance (D_ab = 1 - C_ab/min(N_a,N_b)
    for(my $i=0; $i < @DistMat; $i++){
        for(my $j=$i+1; $j < @DistMat;$j++){
            my ($minCount) = sort {$a <=> $b} @txidCount{@txidList[($i,$j)]};
            my $val = ($minCount > 0) ? $DistMat[$i]->[$j]/$minCount : 0;
            $DistMat[$i]->[$j] = 1 - $val;
            $DistMat[$j]->[$i] = 1 - $val;
        }
    }
    ##Construct 
    my $phylipStr = " ".scalar(@txidList)."\n";
    for (my $i = 0; $i < @txidList; $i++){
        $phylipStr .= sprintf("%-12s%s\n",$txidList[$i],join(" ",@{$DistMat[$i]}));
    }
    my $io = IO::String->new($phylipStr);
    my $matObj = Bio::Matrix::IO->new(-format => "phylip", -fh => $io)->next_matrix;
    my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
    my $dendObj =  $dfactory->make_tree($matObj);
    ReRootTree($dendObj);
    foreach my $node ($dendObj->get_nodes){
        $node->branch_length(undef)
    }
    ##Ensure all nodes in the tree have an id
    my $nextID = 0;
    foreach my $node ($dendObj->get_nodes){
        next if(defined $node->id);
        my $txidKey = join(";",sort (map {$_->id} (get_leaf_descendents($node))));
        if(exists $INodeDict{$txidKey}){
            $node->id($INodeDict{$txidKey});
        } else{
            $node->id(sprintf("IN_%0*d",$_ID_PADDING,$nextID++));
        }
    }
    return $dendObj;
}

#Given a neighbour joining Tree, re-roots the tree on the longest branch
#Input: a Bio::Tree object
#Output: None, Edits the tree Obj in place
sub ReRootTree($){
    my ($treeObj) = @_;
    my $root = $treeObj->get_root_node;
    my @children = $root->each_Descendent;
    return if(scalar(@children) < 2); #Tree is only the root and up to one child, no other way to organize
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

#Returns the makimum number of processors on a linux based system
#Output: -1 if /proc/cpuinfo not found, the number of processors otherwise
sub GetProcessorCount(){
    my $cpu_count = -1;
    if(open my $handle, "/proc/cpuinfo"){
        $cpu_count = scalar(map /^processor/, <$handle>);
        close($handle);
    }
    return $cpu_count;
}

#Given a node in a tree, returns all leaf nodes descendent from that node
#Input: a Bio::Tree::Node objects
#Output: an array of Bio::Tree::Node objects
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


#Utillity function to get the unique elements of an array
#Input: an array
#Output: an array
sub unique(@){
    my %Hash;
    @Hash{@_} = (1) x @_;
    return keys %Hash;
}

#Given a dendrogram and a list of clusters, updates the assignment and penetrance members of the cluster
#Input: a Bio::Tree objects, and an array of CLUSTER objects
#Output: None, CLUSTER objects are modified in place
sub AssignTaxa($$){
    my ($dendObj,$clustListRef) = @_;
    my $pm = Parallel::ForkManager->new($max_proc);
    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $cor_dump, $data_structure_reference) = @_;
            if($exit_code == 0){
                if (defined $data_structure_reference){
                    while( my ($index,$infoRef) = each %{$data_structure_reference}){
                        my($lcaID,$penetrance) = @{$infoRef};
                        $clustListRef->[$index]->assignment($lcaID);
                        $clustListRef->[$index]->penetrance($penetrance);
                    }
                } else {
                    warn "[WARNING] No result from $ident\n";
                }
            } else {
                warn "[WARNING] Unexpected result from block beginning at $ident\n";
            }
        }
    );

    my $nClusters = scalar(@{$clustListRef});
    my $blockSize = ($nClusters < $max_proc) ? 1 : int($nClusters / $max_proc);
    my %lcaDict; # Hash with keys of txidKey, and values of array ref with elements of txid and penetrance
    BLOCK:
    for(my $blockStart = 0; $blockStart < $nClusters; $blockStart += $blockSize){
        $pm->start($blockStart) and next BLOCK;
        my $blockEnd = $blockStart + $blockSize-1;
        $blockEnd = ($blockEnd > ($nClusters-1)) ? ($nClusters - 1) : $blockEnd;
        my %Results;
        for(my $i = $blockStart; $i <= $blockEnd; $i++){
            my $clustObj = $clustListRef->[$i];
            my @txidList = sort unique(@{$clustObj->members});
            my $txidKey = join(";",@txidList);
            my $lcaID = $txidList[0];
            my $penetrance = 100;
            if(exists $lcaDict{$txidKey}){
                ($lcaID,$penetrance) = @{$lcaDict{$txidKey}};
            } else { 
                my @nodeList;
                push(@nodeList,$dendObj->find_node(-id => $_)) foreach (@txidList); 
                my $leafCount = 1;
                if(@nodeList >= 2){
                    my $lcaNode = $dendObj->get_lca(-nodes => \@nodeList);
                    $lcaID = $lcaNode->id;
                    my @leafList = get_leaf_descendents($lcaNode);
                    $leafCount = scalar(@leafList);
                    $penetrance = @nodeList/$leafCount*100
                }
            }
            $lcaDict{$txidKey} = [$lcaID,$penetrance];
            $Results{$i} = [$lcaID,$penetrance];
        }
        $pm->finish(0,\%Results);
    }
    $pm->wait_all_children;
}

            #
###Process Genes
#my $counter = 0;
#my %lcaDict;
#my @geneList = keys(%TxIDSetDict);
#my $geneCount = scalar(keys %TxIDSetDict);
#my $blockSize = int($geneCount / $max_proc);
#$blockSize = 1 if ($blockSize < 1);
#my $targetUpdateTime = 10 * $max_proc;
#GENE:
#for (my $blockStart = 0; $blockStart < @geneList; $blockStart += $blockSize){
#    $pm->start($blockStart) and next GENE;
#    $blockEnd = ($blockEnd > $#geneList) ? $#geneList : $blockEnd;
#    my @genesInBlock = @geneList[$blockStart .. $blockEnd];
#    my @Output;
#    my $lastUpdate = 0;
#    my $count = 0;
#    my $updateRate = 10; #Per thousand progress per update
#    #To introduce an offset for each block to make reporting come consistently rather than in blocks
#    my $initUpdateOffset = $targetUpdateTime - $updateRate * int($blockStart/$blockSize);
#    my $lasttime = time;
#    foreach my $gene (@genesInBlock){
#        $count++;
#        my @txidList = unique(@{$TxIDSetDict{$gene}});
#        my $txidKey = join(';',@txidList);
#        my $lcaID = $txidList[0];
#        my $penetrance = 100;
#        if(exists $lcaDict{$txidKey}){
#            ($lcaID,$penetrance) = @{$lcaDict{$txidKey}};
#        } else {
#            my @nodeList;
#            push(@nodeList,$treeObj->find_node(-id => $_)) foreach (@txidList); 
#            my $leafCount = 1;
#            if(@nodeList >= 2){
#                my $lcaNode = $treeObj->get_lca(-nodes => \@nodeList);
#                $lcaID = $lcaNode->id;
#                my @leafList = get_leaf_descendents($lcaNode);
#                $leafCount = scalar(@leafList);
#                $penetrance = @nodeList/$leafCount*100
#            }
#            $lcaDict{$txidKey} = [$lcaID,$penetrance];
#        }
#        my $progress = int(($count / $blockSize) * 1000);
#        if($bVerbose and $progress - $lastUpdate >= $updateRate){
#            my $curtime = time;
#            my $elapsed = $curtime - $lasttime;
#            $elapsed = 1 if($elapsed < 1);
#            if($elapsed != $targetUpdateTime){
#                $updateRate = ($targetUpdateTime - $initUpdateOffset) * ($progress - $lastUpdate) / $elapsed;
#                $updateRate = 1 if ($updateRate) < 1;
#            }
#            $lastUpdate = $progress;
#            print STDERR "[INFO] The block beginning at $blockStart is @{[$lastUpdate/10]}% complete (@{[$curtime - $lasttime]}s)\n";
#            $lasttime = $curtime;
#            $initUpdateOffset = 0;
#        }
#        my $str = join("\t",($gene,$lcaID,sprintf("%0.2f",$penetrance)));
#        push(@Output,$str);
#    }
#    $pm->finish(0,\@Output);
#}
#$pm->wait_all_children;
#
#
#
#
