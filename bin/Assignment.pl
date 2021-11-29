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
use HUBDesign::Grouping;
use HUBDesign::Logger;

#============================================================
#Declarations

struct  (CLUSTER => {uid => '$', members => '@', assignment => '$', penetrance => '$'});

sub LogParameters();	        #Usage: LogParameters();
sub LoadClusterInfo($);         #Usage: my @clusterList = LoadClusterInfo($clustFile);
sub CountGroupings($@);         #Usage: my @groupingList = CountClusterGroupings($gTreeFile,@cList);
sub EnsureCompatability(\@);    #Usage: my EnsureCompatability(@groupingList);
sub ConstructDendrogram(\@\@);  #Usage: my $treeObj = ConstructDendrogram(@cList,@gist);
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
getopts('g:G:w:t:C:vh',\%opts);
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
my @groupingList = CountGroupings($opts{g},@clusterList);
EnsureCompatability(@groupingList);
my $dendObj = ConstructDendrogram(@clusterList,@groupingList);
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
    my @paramNames = ("Cluster Info","Guide Tree","Dendrogram Output","Cluster Weight","Threads");
    my %params;
    @params{@paramNames} = ($FileDict{ClustInfo},$opts{g},$opts{G},$opts{w},$opts{t});
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

#Given a list of clusters, returns a list of GROUPING objects implied by the list
#Input: an array of CLUSTER Objects
#Output: an array of HUBDesign::Grouping Objects
sub CountGroupings($@){
    $Logger->Log("Counting groupings...","INFO");
    my ($guideTreeFile,@clusterList) = @_;
    my %taxaSet; #Unique Set of taxa in all clusters
    foreach my $clustObj (@clusterList){
        @taxaSet{@{$clustObj->members}} = (1) x scalar(@{$clustObj->members});
    }
    @taxaSet{keys %taxaSet} = (1 .. scalar(keys %taxaSet)); #Assign a unique value to each taxa
    my @groupingList;
    my %GroupingIndex; #Tracks which element of the list a particular grouping is
    foreach my $clustObj (@clusterList){
        #Count the members of the clusters as one group, and all those outside as another
        for (my $inSet = 0; $inSet < 2; $inSet++){
            my @G = unique(@{$clustObj->members});
            @G = notkeys(%taxaSet,@G) unless($inSet);
            @G = sort @G;
            my $id = join(';',sort {$a <=> $b} @taxaSet{@G});
            if(exists $GroupingIndex{$id}){
                $groupingList[$GroupingIndex{$id}]->increment($opts{w});
                $groupingList[$GroupingIndex{$id}]->isObserved(1) if($inSet);
            } else {
                push(@groupingList,HUBDesign::Grouping->new(members => \@G, count => $opts{w},bObs => $inSet));
                $GroupingIndex{$id} = $#groupingList;
            }
        }
    }
    my %SeenSubsets; #Tracks groupings which have already been seen in the tree
    #^ Used to prevent double counting in a cases where non_cluster leaves are in the
    # guide tree
    if(defined $guideTreeFile){
        my $treeObj = Bio::TreeIO->new(-file => $guideTreeFile, -format => "newick")->next_tree();
        my $root = $treeObj->get_root_node;
        my $nextID = 0;
        foreach my $node ($root->get_all_Descendents){
            my @members; 
            foreach my $desc ($node->get_all_Descendents){
                if ($desc->is_Leaf and defined $desc->id and exists $taxaSet{$desc->id}){
                    push(@members,$desc->id);
                }
            }
            next unless (@members);
            my $id = join(';',sort {$a <=> $b} @taxaSet{@members});
            next if(exists $SeenSubsets{$id});
            $SeenSubsets{$id}=1;
            if(exists $GroupingIndex{$id}){
                $groupingList[$GroupingIndex{$id}]->increment;
                $groupingList[$GroupingIndex{$id}]->isObserved(1);
            } else {
                push(@groupingList,HUBDesign::Grouping->new(members => \@members, count => 1,bObs=>1));
                $GroupingIndex{$id} = $#groupingList;
            }
        }
    }
    $Logger->Log(sprintf("Counted occurences for %d groupings",scalar(@groupingList)),"INFO");
    return @groupingList;
}



#Given an ordered list of groupings, remove any which are incompatible with earlier groupings in the
#tree #smaller is a subset of the larger)
#Input:  Reference to an array of ordered HUBDesign::Grouping objects
#Output: None, edits the input in place
sub EnsureCompatability(\@){
    $Logger->Log("Culling incompatible groupings...","INFO");
    my $gListRef = shift;
    #Sort Groupings By Number of observances,
    #Then prefer groupings actually observed, rather than implied inverting a cluster
    #Then break ties by the size of the grouping
    @{$gListRef} = sort {$b->count <=> $a->count || $b->isObserved <=> $a->isObserved || scalar($b->members) <=> scalar($a->members)} @{$gListRef};
    for(my $i = 0; $i < $#{$gListRef}; $i++){
        for(my $j = $i+1; $j < @{$gListRef}; $j++){
            my $bCompat = $gListRef->[$i]->isCompatible($gListRef->[$j]);
            #print STDERR join("-",$gListRef->[$i]->members),"\t",join("-",$gListRef->[$j]->members),"\t$bCompat\n";
            next if($bCompat);
            splice(@{$gListRef},$j--,1);
        }
    }
    $Logger->Log(sprintf("A total of %d compatible groupings remain",scalar(@{$gListRef})),"INFO");
}


#Given a set of clusters, and a set of compatible groupings of those clusters, generates a
#tree object 
#Input: a set of clusters (used to ensure all input taxa are in the tree), and a set of
#compatible groupings
#Output: a Bio::Tree Object
sub ConstructDendrogram(\@\@){
    $Logger->Log("Constructing Dendrogram...","INFO");
    my ($cListRef,$gListRef) = @_;
    #Sort the groupings by the coverage of the set
    @{$gListRef} = sort {scalar($b->members) <=> scalar($a->members)} @{$gListRef};
    my $firstSingleton = 0;
    while($firstSingleton < @{$gListRef} and scalar($gListRef->[$firstSingleton]->members) > 1){
        $firstSingleton++ 
    }
    splice(@{$gListRef},$firstSingleton) if($firstSingleton < @{$gListRef});
    #Array of array ref with index elements; tracks which groupings are directly nested inside another; Note this array has one extra element; the first element is the children of the root
    my @children = ([]);
    push(@children,[]) foreach (@{$gListRef});
    for(my $i = 0; $i < @{$gListRef}; $i++){
        my $j = 0;
        my $bSubset = 1;
        #print STDERR "$i\n";
        while($bSubset){ #Until the current grouping is not a subset of any other child of the current node
            #  print "#$i\t$j\t$bSubset\n";
            $bSubset=0;
            #Find which child (if any) in the current node of which the current group is a subset
            foreach my $cIdx (@{$children[$j]}){
                if($gListRef->[$i]->isSubset($gListRef->[$cIdx])){
                    $bSubset=1;
                    $j = $cIdx+1;
                    last;
                }
            }
        }
        push(@{$children[$j]},$i);
    }
    my %taxaSet; #Unique Set of taxa which must eventually have a leaf in the current node
    #Initialize this with all taxa at the root
    foreach my $clustObj (@{$cListRef}){
        @taxaSet{@{$clustObj->members}} = (1) x scalar(@{$clustObj->members});
    }
    my @nodeList;
    push(@nodeList,Bio::Tree::Node->new()) foreach (@children);
    for(my $nIdx = 0; $nIdx < @children; $nIdx++){
        my $node = $nodeList[$nIdx];
        unless($nIdx == 0){
            %taxaSet = ();
            @taxaSet{$gListRef->[$nIdx-1]->members} = (1) x scalar($gListRef->[$nIdx-1]->members);
        }
        #Add any leaves which should be within the tree under this node, but are not in
        #any of the groupings marked as children for this node
        foreach my $txid (notkeys(%taxaSet,map {$gListRef->[$_]->members} @{$children[$nIdx]})){
            $node->add_Descendent(Bio::Tree::Node->new(-id => $txid));
        }
        foreach my $cIdx (@{$children[$nIdx]}){
            if(scalar($gListRef->[$cIdx]->members) == 1){
                $nodeList[$cIdx+1]->id(($gListRef->[$cIdx]->members)[0]);
            }
            $node->add_Descendent($nodeList[$cIdx+1]);
        }
    }
    my $dendObj = Bio::Tree::Tree->new(-root => $nodeList[0]);
    ##Ensure all nodes in the tree have an id
    my $nextID = 0;
    foreach my $node ($dendObj->get_nodes){
        next if(defined $node->id);
        $node->id(sprintf("IN_%0*d",$_ID_PADDING,$nextID++));
    }
    Bio::TreeIO->new(-format => "newick",-file => ">".$opts{G})->write_tree($dendObj);
    $Logger->Log("Dendrogram written to $opts{G}","INFO");
    return $dendObj;
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
