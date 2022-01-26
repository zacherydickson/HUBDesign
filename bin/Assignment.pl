#!/usr/bin/env perl
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
use FindBin;
use lib File::Spec->catdir($FindBin::RealBin, '..', 'lib');
use HUBDesign::Util qw(ValidateThreadCount ProcessNumericOption OpenFileHandle LoadConfig);
use HUBDesign::Bipartition;
use HUBDesign::Logger;

#============================================================
#Declarations

struct  (CLUSTER => {uid => '$', members => '@', assignment => '$', penetrance => '$'});

sub LogParameters();	        #Usage: LogParameters();
sub LoadClusterInfo($);         #Usage: my @clusterList = LoadClusterInfo($clustFile);
sub FastGetLeaves ($$);         #Usage: FastGetLeaves($node,$leaveListRef);
sub AddBip($$$$@);		#Usage: AddBip(\%IndexMap,\%excludeSet,\%taxaMap,$weight,@taxaList);
sub CountBips($@);         	#Usage: CountBips($guideTreeFile,@clusterList);
sub EnsureCompatability(\@);    #Usage: my EnsureCompatability(@groupingList);
sub ProcessNode($$$$@);		#Usage: ProcessNode($node,$gListRef,$indexMapRef,$bFirst,@indexList);
sub ConstructDendrogram(\@\@@);	#Usage: ConstructDendrogram($cListRef,$gListRef,@outgrp);
sub GetProcessorCount();        #Usage: my $sys_proc_count = GetProcessorCount();
sub get_leaf_descendents($);    #Usage: my @leafList = get_leaf_descendents($lcaNode);
sub unique(@);                  #Usage: my @unique = unique(@list)
sub notkeys(\%@);               #Usage: my @missedKeys = notkeys(%hash,@keys);
sub AssignTaxa($$);             #Usage: AssignTaxa($dendObj,\@clusterList);
sub ParseConfig($); 	        #Usage: ParseConfig($ConfigFile);

my $_ID_PADDING = 5;
my %DEFAULT = (G => 'Guide.dend', w => 1, t => 1);

#============================================================
#Handling Input

my %opts;
getopts('g:G:w:t:C:Ovh',\%opts);
ParseConfig($opts{C});

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
            "-G PATH\tThe file to which to write the dendrogram used to assign clusters to taxa (Default: $DEFAULT{G})\n".
            "-w (0,∞]\tThe relative weight of the topology of the Clustering compared to the\n".
            "\t\tGuide Tree (Default: $DEFAULT{w})\n".

            "-t INT\tNumber of Threads; 0 = sys_max (Default: $DEFAULT{t})\n".
            "-C PATH\t Path to a HUBDesign Config file (Default: $FindBin::RealBin/../HUBDesign.cfg)\n".
            "===Flags\n".
            "-O\tDisable outgrouping on most isolated taxa\n".
            "-v\tVerbose Output\n".
            "-h\tDisplay this message and exit\n";
    } else{
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

foreach my $opt (keys %DEFAULT){
    $opts{$opt} = $DEFAULT{$opt} unless(exists $opts{$opt});
}

#Init Logger Verbosity
$opts{v} = (exists $opts{v}) ? "INFO" : "WARNING";
my $Logger = HUBDesign::Logger->new(level => $opts{v});

#Process Simple Numeric Options
$opts{t} = ProcessNumericOption($opts{t},$DEFAULT{t},0,undef,1,"Max Threads");
$opts{w} = ProcessNumericOption($opts{w},$DEFAULT{w},0,undef,0,"Guide Tree Weight");

if(exists $opts{g}){
    $Logger->Log("Could not locate file: $opts{g}",'ERROR') unless(-e $opts{g});
} else{
    $opts{g} = undef;
}
$opts{G} = (exists $opts{G}) ? $opts{G} : "Guide.dend";
my $max_proc = ValidateThreadCount($opts{t});

my %FileDict;
@FileDict{qw(ClustInfo)} = @ARGV[0 .. 0];
while( my ($type,$file) = each %FileDict){
    $Logger->Log("$type file ($file) does not exist","ERROR") unless(-e $file);
}

#============================================================
#Main Script

LogParameters();
my @clusterList = LoadClusterInfo($FileDict{ClustInfo});
my @groupingList;
my @outgroupTaxa = CountBips($opts{g},@clusterList);
EnsureCompatability(@groupingList);
my $dendObj = ConstructDendrogram(@clusterList,@groupingList,@outgroupTaxa);
AssignTaxa($dendObj,\@clusterList);
print join("\t",qw(ClusterID TaxonID Penetrance)),"\n";
foreach my $clustObj (@clusterList){
    print join("\t",($clustObj->uid,$clustObj->assignment,sprintf("%.2f",$clustObj->penetrance))),"\n";
}
$Logger->Log("Done","INFO");


#============================================================
#Subroutines

#Outpits the Parameters to STDERR
sub LogParameters(){
    my @paramNames = ("Cluster Info","Guide Tree","Dendrogram Output","Cluster Weight","Outgrouping","Threads");
    my %params;
    @params{@paramNames} = ($FileDict{ClustInfo},$opts{g},$opts{G},$opts{w},"Enabled",$opts{t});
    $params{Outgrouping} = 'Disabled' if(exists $opts{O});
    delete $params{"Guide Tree"} unless(defined $opts{g});
    $Logger->LogParameters(map {($_,$params{$_})} @paramNames);
}

#Given a cluster info file, parse out the cluster ids and their associated taxon ids
#Input: a path to a cluster info file
#Output: an array of CLUSTER objects
sub LoadClusterInfo($){
    $Logger->Log("Loading Cluster Info...","INFO");
    my $fh = OpenFileHandle(shift,"ClustInfo","ERROR");
    my @clusterList;
    my %indexDict;
    my @headers = split(/\t/,<$fh>);
    if(@headers < 2){
        $Logger->Log("ClustInfo File must have 2 tab delimited columns","ERROR");
    }
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
    $Logger->Log(sprintf("Loaded %d clusters",scalar(@clusterList)),"INFO");
    return(@clusterList);
}

#Depth First traversal of a tree getting all tips decendent from each node
# Input: A node and a reference to an array to store the results in
# Output: The list of tip ids descendent from the node
sub FastGetLeaves ($$){
    my $node = shift;
    my $outputRef = shift;
    if($node->is_Leaf){
        push(@{$outputRef},[$node->id]);
        return ($node->id);
    }
    my @taxaList; 
    push(@{$outputRef},[]);
    my $index = $#{$outputRef};
    foreach my $child ($node->each_Descendent){
        my @ctaxaList = FastGetLeaves($child,$outputRef);
        push(@{$outputRef->[$index]},@ctaxaList);
    }
    return @{$outputRef->[$index]};
}


#Input: A reference to an array of groupings Objects
#	A reference to a hash with keys of ids and values of indexes in the first argument
#	A reference to a hash with keys of grouping ids which shouldn't be included
#	A reference to a hash with keys of taxids and values of indexes for the taxids
#	An array of taxids
#Ouput: A HUBDesign::Bipartition Object constructed from the taxaList
sub AddBip($$$$@) {
    my ($gIndexRef,$exSetRef,$txMapRef,$weight,@taxaList) = @_;
    my @indexList = @{$txMapRef}{@taxaList};
    my $grpObj = HUBDesign::Bipartition->new(members => \@indexList, count => $weight);
    my $id = $grpObj->get_bits;
    return if(exists $exSetRef->{$id});
    $exSetRef->{$id} = 1;
    if(exists $gIndexRef->{$id}){
        $groupingList[$gIndexRef->{$id}]->increment($weight);
    } else {
        push(@groupingList,$grpObj);
        $gIndexRef->{$id} = $#groupingList;
    }
    return $grpObj;
}

#Given a list of clusters, returns a list of GROUPING objects implied by the list
#Input: an array of CLUSTER Objects
#Output: an array of outgroup taxa
sub CountBips($@){
    my ($guideTreeFile,@clusterList) = @_;
    $Logger->Log("Counting Bipartition Instances...","INFO");
    my %taxaSet; #Unique Set of taxa in all clusters
    foreach my $clustObj (@clusterList){
        @taxaSet{@{$clustObj->members}} = (1) x scalar(@{$clustObj->members});
    }
    my $nTaxa = scalar(keys %taxaSet);
    HUBDesign::Bipartition->initialize_nTaxa($nTaxa);
    @taxaSet{sort keys %taxaSet} = (0 .. ($nTaxa-1)); #Assign a unique value to each taxa
    my %taxaCount; #For tracking how many times each taxa actually appears in any non singular
    @taxaCount{keys %taxaSet} = (0) x $nTaxa;
    my %GroupingIndex; #Tracks which element of the list a particular grouping is
    foreach my $clustObj (@clusterList){
        my @taxaList = unique(@{$clustObj->members});
        if(@taxaList > 1){
            foreach my $txid (@taxaList){
                $taxaCount{$txid}++;
            }
        }
        my $grpObj = AddBip(\%GroupingIndex,{},\%taxaSet,$opts{w},@taxaList);
    }
    my %SeenSubsets; #Tracks groupings which have already been seen in the tree
    #^ Used to prevent double counting in cases where non-cluster leaves are in the
    # guide tree
    if(defined $guideTreeFile){
        my $treeObj = Bio::TreeIO->new(-file => $guideTreeFile, -format => "newick")->next_tree();
        my $root = $treeObj->get_root_node;
        my @leavesDescendentFromEachNode;
        FastGetLeaves($root,\@leavesDescendentFromEachNode);
        shift @leavesDescendentFromEachNode; #Don't want to count the bipartition at the root of All|Empty
        if(scalar($root->each_Descendent) == 2){ #Two children at the root are complementary, only one is needed
            shift @leavesDescendentFromEachNode;
        }
        foreach my $taxaListRef (@leavesDescendentFromEachNode){
            my @taxaList;
            foreach my $txid (@{$taxaListRef}){
                push(@taxaList,$txid) if(exists $taxaSet{$txid});
            }
            next unless(@taxaList);
            if(@taxaList > 1){
                foreach my $txid (@taxaList){
                    $taxaCount{$txid}++;
                }
                AddBip(\%GroupingIndex,\%SeenSubsets,\%taxaSet,1,@taxaList);
            }
        }
    }
    my $bipCount = scalar(@groupingList);
    #Find the Taxa which appears in the fewest groupings
    my %Min = (val => $bipCount, txidList => []);
    while(my ($txid,$count) = each %taxaCount){
        if($count < $Min{val}){
            %Min = (val => $count, txidList => [$txid]);
        } elsif($count == $Min{val}){
            push(@{$Min{txidList}}, $txid);
        }
    }
    $Logger->Log("Instances counted for $bipCount Bipartitions","INFO");
    return @{$Min{txidList}};
}

#Given an ordered list of groupings, remove any which are incompatible with earlier groupings in the tree
#Input:  Reference to an array of HUBDesign::Bipartition objects
#Output: None, edits the input in place
sub EnsureCompatability(\@){
    $Logger->Log("Culling Incompatibilities...","INFO");
    my $gListRef = shift;
    my $nBip = scalar(@{$gListRef});
    #Sort Groupings By Number of observances, then break ties by the size of the grouping
    @{$gListRef} = sort {$b->count <=> $a->count || $b->size <=> $a->size} @{$gListRef};
    my $comparisons = $nBip * ($nBip - 1) / 2;
    my $elim = 0;
    for(my $i = 0; $i < $#{$gListRef}; $i++){
        for(my $j = $i+1; $j < @{$gListRef}; $j++){
            my $bCompat = $gListRef->[$i]->isCompatible($gListRef->[$j]);
            next if($bCompat);
            splice(@{$gListRef},$j--,1);
        }
    }
    $nBip = scalar(@{$gListRef});
    $Logger->Log("There are $nBip compatible bipartitions remaining ","INFO");
}

# Given a set of indexes to be covered by a node, both recursively and
# iteratively adds nodes until all the required children are decendents of the
# node
# Input:    A Bio::Tree::Node object describing the current node
#           A reference to an array of HUBDesign::Bipartition Objects
#           A reference to a hash with keys of indexes and values of taxon ids
#           A list of indexes to be covered by the node
# Ouput:    None, alll nodes are attached to the node node objects 
sub ProcessNode($$$$@){
    my ($node,$gListRef,$indexMapRef,$bFirst,@nodeIndexList) = @_;
    @nodeIndexList = sort {$a <=> $b} @nodeIndexList;
    my $nTaxa = scalar(keys %{$indexMapRef});
    while (@nodeIndexList){
        my $nodeGrpSize = scalar(@nodeIndexList);
        if($nodeGrpSize == 1){ #Is tip
            $node->id($indexMapRef->{shift @nodeIndexList});
            return;
        }
        my $flag = ($nodeGrpSize > $nTaxa/2) ? 2 : 0; #Make sure comparisons are always to the 'on' grp
        #Find the bip which contains a grouping which is the largest subset of the taxaSet
        my %Max = (size => -1,obj => undef,isNeg => -1,idx => -1);
        my $nodeGrp = HUBDesign::Bipartition->new(members => \@nodeIndexList);
        if($bFirst){
            %Max = (size => $gListRef->[0]->size,obj => $gListRef->[0],isNeg => 1, idx => 0);
            $bFirst = 0;
        } else {
            for(my $i = 0; $i < @{$gListRef};$i++){
                my $grpObj = $gListRef->[$i];
                next if($grpObj->size > $nodeGrpSize);
                my %cur; 
                #Compare Big grp in bip to nodeGrp
                if($nTaxa-$grpObj->size < $nodeGrpSize && $grpObj->isSubset($nodeGrp,$flag | 0b1)){
                    %cur = (size => $nTaxa-$grpObj->size,obj => $grpObj, isNeg => 1, idx => $i);
                }
                #Compare Small grp in bip to nodeGrp
                if(!%cur && $grpObj->size < $nodeGrpSize && $grpObj->isSubset($nodeGrp,$flag | 0b0)){
                    %cur = (size => $grpObj->size, obj => $grpObj, isNeg => 0, idx => $i);
                }
                %Max = %cur if(%cur && $cur{size} > $Max{size});
            }
        }
        if(defined $Max{obj}){ #Found a subset
            #Remove the found bip from further consideration
            splice(@{$gListRef},$Max{idx},1);
            #Make new nodes for the subset and the remainder
            my @children;
            foreach my $i (0 .. 1){
                push(@children,Bio::Tree::Node->new());
                $node->add_Descendent($children[$i]);
            }
            $node = $children[1];
            #Recurse with the subset
            my @indexList = $Max{obj}->members;
            @indexList = notkeys(%{$indexMapRef},@indexList) if ($Max{isNeg});
            ProcessNode($children[0],$gListRef,$indexMapRef,0,@indexList);
            #Repeat the loop with the remainder
            my %indexSet;
            @indexSet{@nodeIndexList} = (1) x $nodeGrpSize;
            @nodeIndexList = notkeys(%indexSet,@indexList);
            unless(@nodeIndexList){
                $children[0]->ancestor($node->ancestor->ancestor);
                $node->ancestor->remove_Descendent($node);
            }
        } else{ #This is an internal node, but there are no subgroupings, add children for each remaining taxa
            foreach my $index (@nodeIndexList){
                $node->add_Descendent(Bio::Tree::Node->new(-id => $indexMapRef->{$index}))
            }
            @nodeIndexList = ();
        }
    }
}

#Given a set of compatible bipartitions of a tree, constructs the tree
#Works by recursively subdividing the entire set of taxa
#Input: A reference to an array of CLUSTER objects
#       A reference to an array of HUBDesign::Bipartition Objects
#       A list of taxa which are to be considered outgroups
#Output:    A Bio::Tree Object
sub ConstructDendrogram(\@\@@){
    $Logger->Log("Constructing Dendrogram...","INFO");
    my ($cListRef,$gListRef,@outgroupTaxa) = @_;
    #Construct a mapping from bitstr index to txid
    my %taxaSet;
    foreach my $clustObj (@{$cListRef}){
        @taxaSet{@{$clustObj->members}} = (1) x scalar(@{$clustObj->members});
    }
    my $nTaxa = scalar(keys %taxaSet);
    my %indexMap;
    @indexMap{(0 .. ($nTaxa-1))} = sort keys %taxaSet;
    my $root = Bio::Tree::Node->new();
    ProcessNode($root,$gListRef,\%indexMap,1,(keys %indexMap));
    my $dendObj = Bio::Tree::Tree->new(-root => $root);
    if(!exists $opts{O} and @outgroupTaxa){
        my @outgroupNodes;
        foreach my $txid (@outgroupTaxa){
            push(@outgroupNodes,$dendObj->find_node(-id => $txid));
        }
        $Logger->Log("Rooting with outgroup(s): @outgroupTaxa...","INFO");
        if(@outgroupNodes){
            my $outgroupNode = $outgroupNodes[0];
            $outgroupNode = $dendObj->get_lca(@outgroupNodes) if(@outgroupNodes > 1);
            if(defined $outgroupNode->ancestor and defined $outgroupNode->ancestor->ancestor){
                $dendObj->reroot($outgroupNode->ancestor);
                $dendObj->splice($root);
                foreach my $messNode ($dendObj->find_node(-branch_length => 0)){
                    $dendObj->splice($messNode);
                }
            } else {
                $Logger->Log("New root is the old root","WARNING");
            }
        } else {
            $Logger->Log("Could not find @outgroupTaxa in the tree : Skipping Outgrouping","WARNING");
        }
    }
    #Ensure Every Node in the Tree Has an ID
    my $next_id = 0;
    foreach my $node ($dendObj->get_nodes){
        $node->id(sprintf("IN_%0*d",$_ID_PADDING,$next_id++)) unless(defined $node->id);
    }
    my $treeIO = Bio::TreeIO->new(-file => ">$opts{G}", -format => 'newick');
    $treeIO->write_tree($dendObj);
    $Logger->Log("Dendrogram written to $opts{G}","INFO");
    return $dendObj;
}

#Returns the maximum number of processors on a linux based system
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

#Utility  function to get all keys of a hash except for the provided set
#Input a hash reference, and an array of keys
#Output an array
sub notkeys(\%@){
    my ($hashref,@keys) = @_;
    my %tmpHash = %$hashref;
    delete @tmpHash{@keys};
    keys %tmpHash;
}

#Given a dendrogram and a list of clusters, updates the assignment and penetrance members of the cluster
#Input: a Bio::Tree objects, and an array of CLUSTER objects
#Output: None, CLUSTER objects are modified in place
sub AssignTaxa($$){
    $Logger->Log("Assigning Taxa...","INFO");
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
            my @txidList = sort (unique(@{$clustObj->members}));
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
                    if($lcaNode->is_Leaf){
                        warn "[WARNING] LCA ($lcaID) for ".scalar(@nodeList)." nodes (".join(",",@txidList).") found to be a leaf\n";
                    } else {
                        my @leafList = get_leaf_descendents($lcaNode);
                        $leafCount = scalar(@leafList);
                    }
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

sub ParseConfig($){
    my $file = shift;
    $file = "$FindBin::RealBin/../HUBDesign.cfg" unless defined $file;
    unless(-e $file){
        warn "Could not find Config file ($file): Revert to Hardcoded defaults\n";
        return;
    }
    my %Config = LoadConfig($file,"WARNING");
    $_ID_PADDING = $Config{'ID-padding'} if(exists $Config{'ID-padding'});
    $DEFAULT{t} = $Config{'threads'} if(exists $Config{'threads'});
    $DEFAULT{w} = $Config{'weight-guide-tree'} if(exists $Config{'weight-guide-tree'});
}
