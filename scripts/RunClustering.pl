#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;
use Getopt::Std;
use Class::Struct;
use File::Temp;

#============================================================
#Declations

struct  (GENE => {tid => '$', chr => '$', uid => '$', start => '$', end => '$', strand => '$', seq => '$'});
struct  (CLUSTER => {uid => '$', seq => '$', members => '@'});

my $CLUST_ID_PADDING = 5;
my $TEMP_DIR = "/tmp";
my $ALIGNER = "mafft --auto --localpair --quiet";
my $ALIGNER_FORMAT = "fasta";

sub LoadGFFSet($@);         #Usage: my %GFFSet = LoadGFFSet($GFFListFile,@GFFFileList);
sub LoadGeneInfo(%);        #Usage: my %GeneInfo = LoadGenes(%GFFSet);
sub HandleSinglets($$$);    #Usage: my $nSinglets = HandleSinglets(\%GeneInfo,$ClustHandle,$SeqOObj);
sub OutputCluster($$$);     #Usage: OutputCluster($ClustObj,$ClustHandle,$SeqOObj);
sub AlignByProt($$);        #Usage: AlignByProt(\@geneList,$geneName);
sub ConstructTree($$);      #Usage: my $TreeObj = ConstructTree($alnObj,$geneName);
sub ReRootTree($);          #Usage: ReRootTree($treeObj);
sub ExtractClusters($$$$);  #Usage: my @clusterList = ExtractClusters(\@geneList,$alnObj,$treeObj,$geneName);
sub PartitionTipLabels($$); #Usage: my @ClusteredLAbels = PartitionTipLAbels($nodeObj,$threshold)

#============================================================
#Handling Input

my %opts;
getopts('d:e:l:p:o:t:vVh',\%opts);

if(@ARGV < 1 or exists $opts{h}){
    my $Usage = basename($0)." [-options] [-l GFF.tsv] GFF1.gff[.gz] ... > ClustSeq.fna";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of GFF files, extracts all CDS's and moves them to a file with the gene name\n".
            "===Usage\n\t$Usage\n".
            "===Input\n".
            "\tGFF1 ...\tPaths to gff files each of which include the sequence at the end\n".
            "\t\tmarked by ##FASTA\n".
            "\t\tIdeally, these GFFs are output from Prokka, and have ID=XXX;Name=XXX in\n".
            "\t\tThe ninth column\n".
            "\t\tNote: The basenames of the gff files are used as taxon IDs\n".
            "\t\t\tIf this is undesirable, use the -l option to specifiy input files\n".
            "===Output\n".
            "\tA fasta formatted file with Cluster IDs has headers\n".
            "===Options\n".
            "-d [0-1]\tThe Maximum root-tip divergence within clusters (Default: 0.15)\n".
            "-e STR\tA comma separated list of file extensions [Default: '.gz,.gff']\n".
            "-l PATH\tA tab delimited file with a path to a gff file, and\n".
            "\t\tan optional taxonID on each line\n".
            "\t\tNote: If -l is specified, then no command line GFF's need to be specified\n".
            "-p INT\tTranslation Table to use (Default: 1 [Standard])\n".
            "-o PATH\t The file to which cluster info will be written (Default: ClustInfo.tsv)\n".
            "-t STR\tA string acting as a flag for the external aligner for the number of threads [Default '']\n".
            "===Flags\n".
            "-v\tVerbose Output\n".
            "-h\tDisplay this message an exits\n";
    } elsif(! exists $opts{l}) {
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}


$opts{d} = (exists $opts{d}) ? $opts{d} : 0.15;
if($opts{d} < 0 or $opts{d} > 1){
    die "Clustering threshold must be a floating value between 0 and 1 inclusive\n";
}
$opts{e} = exists $opts{e} ? [split(",",$opts{e})] : [".gz",".gff"];
if(exists $opts{l}){
    die "Could not locate file: $opts{l}\n" unless(-e $opts{l});
} else{
    $opts{l} = undef;
}
$opts{p} = exists $opts{p} ? $opts{p} : 1;
$opts{o} = exists $opts{o} ? $opts{o} : "ClustInfo.tsv";
$opts{t} = exists $opts{t} ? $opts{t} : "";
$opts{v} = 1 if (exists $opts{V});


#============================================================
#Main Script

my $SeqOObj = Bio::SeqIO->newFh(-format => 'fasta');
open (my $clustHandle, ">$opts{o}") or die "Could not write to $opts{o}: $!\n";
print $clustHandle join("\t",qw(ClusterID TaxonID SeqID GeneID Start End Strand)),"\n";


my %GeneInfo = LoadGeneInfo(LoadGFFSet($opts{l},@ARGV));
printf STDERR "Gene families Loaded:\t%d\n",scalar(keys %GeneInfo) if(exists $opts{v});
my $nSinglets = HandleSinglets(\%GeneInfo,$clustHandle,$SeqOObj);
printf STDERR "Single member gene families:\t%d\n",$nSinglets if(exists $opts{v});
foreach my $geneName (keys %GeneInfo){
    printf STDERR "Processing gene family %s...\n",$geneName if(exists $opts{v});
    my $seqListRef = $GeneInfo{$geneName};
    print STDERR "\tAligning...\n" if(exists $opts{v});
    my $alnObj = AlignByProt($seqListRef,$geneName);
    print STDERR "\tGenerating Tree...\n" if(exists $opts{v});
    my $treeObj = ConstructTree($alnObj,$geneName);
    print STDERR "\tExtracting Clusters...\n" if(exists $opts{v});
    my @clusterList = ExtractClusters($seqListRef,$alnObj,$treeObj,$geneName);
    foreach my $clusterObj (@clusterList){
        OutputCluster($clusterObj,$clustHandle,$SeqOObj);
    }
}
print STDERR "Done\n" if(exists $opts{v});

close($clustHandle);

#============================================================
#Functions

#Handles combinations of gff files provided on the command line,
# and in a list file, with or without taxon ids, and with various file extensions.
#Outputs a hash with keys of taxon IDs and values of gffFilePaths
# %gffSet = (TID => Path, ...)
sub LoadGFFSet($@){
    my $listFile = shift;
    #array with elements of [$file,$tid]
    my @fileList = map {[$_,undef]} @_;
    my %gffSet;
    ##Process GFF List File if provided
    if(defined $listFile){
        open(my $fh, $listFile) or die "Could not open $listFile: $!\n";
        my @lines = <$fh>;
        close($fh);
        chomp(@lines);
        push(@fileList,map {[split(/\t/,$_)]} @lines);
    }
    ##Process Files given on the command line
    while (@fileList){
        my ($file,$tid) = @{shift(@fileList)};
        $tid = basename($file,@{$opts{e}}) unless(defined $tid);
        if(exists $gffSet{$tid}){
            warn "[WARNING] Multiple files with the same taxonID ($tid) detected: ".
                "Ignoring all but $gffSet{$tid}\n";
        } else {
            $gffSet{$tid} = $file;
        }
    }
    return %gffSet;
}

#Parses GFF Files and extracts the location, identifiers, and sequences for every gene in
# the input set
#Takes a hash with keys of taxonIDs and values of gffFilePaths as input i.e Output of
# LoadGFFSet)
#Outputs a hash with keys of Gene Name and values of array ref with GENE elements
# %GeneInfo = (Name => [GENE, ...], ...);
sub LoadGeneInfo(%){
    my %gffSet = @_;
    my %GeneInfo;
    ##Process Each GFF File
    while( my ($tid,$filename) = each %gffSet){
        #Handle potentially gzipped files
        my $opencmd = $filename;
        if ($filename =~ m/\.gz$/){
            $opencmd = "gzip -dc $filename |"
        }
        if(open(my $fh, $opencmd)){
            ##Parse the portion of the GFF with gene info
            my %chrIndex; #For rapidly filling in sequence info
            while(my $line = <$fh>){
                last if($line =~ m/##FASTA/); #Reached Sequence info
                next if($line =~ m/^##/); # Skip headers
                chomp($line);
                my ($seqID,undef, $type, $s, $e, undef, $strand, undef, $infoStr) =
                    split(/\t/,$line);
                next unless($type eq "CDS");
                my %info = split(/=|;/,$infoStr);
                my $geneName = exists $info{Name} ? $info{Name} : "UNKNOWN";
                ($geneName) = split(/_/,$geneName);
                $GeneInfo{$geneName} = [] unless(exists $GeneInfo{$geneName});
                push(@{$GeneInfo{$geneName}},GENE->new(uid => $info{ID},start => $s,
                        end => $e, strand => $strand, seq => undef, tid => $tid, chr => $seqID));
                my $label = join(';',($seqID,$geneName));
                $chrIndex{$label} = [] unless(exists $chrIndex{$label});
                push(@{$chrIndex{$label}},$#{$GeneInfo{$geneName}});
            }
            ##Extract the Sequence Info
            my $SeqI = new Bio::SeqIO(-fh => $fh, -format => 'fasta');
            while(my $seqObj = $SeqI->next_seq){
                my $id = $seqObj->id;
                foreach my $geneName (keys %GeneInfo){
                    my $label = join(';',($id,$geneName));
                    next unless(exists $chrIndex{$label});
                    foreach my $index (@{$chrIndex{$label}}){
                        my $geneObj = $GeneInfo{$geneName}->[$index];
                        my $seq = $seqObj->subseq($geneObj->start,$geneObj->end);
                        if($geneObj->strand eq '-'){
                            $seq = Bio::Seq->new(-seq => $seq, -id => $seqObj->id)->revcom()->seq;
                        }
                        $geneObj->seq($seq);
                    }
                }
            }
            close($fh);
        } else {
            warn "Could not open $filename: $!\n";
        }
    }
    return %GeneInfo;
}

#Outputs any genes which have only one sequence as their own complete cluster
#Takes a reference to a hash with Keys of Gene Names and values of array ref with GENE elements (Output of
#   LoadGeneInfo)
#   Also Takes a file handle and a Bio::SeqIO Obj specifying where to write the info and sequence to
#Returns The number of singlet Genes removed
sub HandleSinglets($$$){
    my ($GeneInfo,$clustHandle,$SeqOObj) = @_;
    my $nSinglets = 0;
    my @geneNameList = keys %{$GeneInfo};
    foreach my $geneName (@geneNameList){
        my $geneListRef = $GeneInfo->{$geneName};
        if(scalar(@{$geneListRef}) == 1){
            my $clustID = sprintf("%s_%0*d",$geneName,$CLUST_ID_PADDING,0);
            my $clustObj = CLUSTER->new(uid => $clustID, seq => $geneListRef->[0]->seq, members => $geneListRef);
            OutputCluster($clustObj,$clustHandle,$SeqOObj);
            delete $GeneInfo->{$geneName};
            $nSinglets++;
        }
    }
    return $nSinglets;
}

#Outputs the info and sequence of a cluster to their appropriate files
#Take a Cluster Obj, and handles for the clust info and sequence files
#Returns Nothing
sub OutputCluster($$$){
    my ($clustObj,$clustHandle,$SeqOObj) = @_;
    my $oldFH = select($clustHandle);
    foreach my $geneObj (@{$clustObj->members}){
        print join("\t",($clustObj->uid,$geneObj->tid,$geneObj->chr,$geneObj->uid,$geneObj->start,$geneObj->end,$geneObj->strand)),"\n";
        print $SeqOObj Bio::Seq->new(-seq => $clustObj->seq, -id => $clustObj->uid);
    }
    select($oldFH);
}

#Aligns a set of sequences by their protein sequence
#Input: a reference to a list of GENE Objects
#Output: a Bio::Align::AlignI compliant object
sub AlignByProt($$){
    my ($geneListRef,$geneName) = @_;
    my @Alignment;
    ## Output Translated Sequences to a temporary File
    my $tmpFile = File::Temp->new(DIR => $TEMP_DIR, TEMPLATE => "HUBCLUSTXXXX", SUFFIX =>".fna",UNLINK=>1);
    print STDERR "\t\tOutputing translated sequences to $tmpFile...\n" if(exists $opts{V});
    my $SeqO = Bio::SeqIO->newFh(-format => "fasta", -fh => $tmpFile);
    my %geneIndex;
    my $nextIndex=0;
    foreach my $geneObj (@{$geneListRef}){
        print $SeqO Bio::Seq->new(-seq => $geneObj->seq, -id => $geneObj->uid)->translate(-codontable_id => $opts{p}, -complete => 1);
        $geneIndex{$geneObj->uid} = $nextIndex++;
    }
    ## Call the external aligner, and get its output
    print STDERR "\t\tAligning translated sequences...\n" if(exists $opts{V});
    open( my $fh, "-|", "$ALIGNER $opts{t} ".$tmpFile->filename) or warn "Alignment Error for $geneName: $!";
    my $AlnI = Bio::AlignIO->new(-fh => $fh, -format => $ALIGNER_FORMAT);
    my $alnObj = $AlnI->next_aln();
    ## Detranslate the protein sequences
    print STDERR "\t\tReturning aligned proteins to nucleotide sequence...\n" if(exists $opts{V});
    foreach my $seqObj ($alnObj->each_seq){
        my $id = $seqObj->id;
        next unless(exists $geneIndex{$id});
        my @aaList = split("",$seqObj->seq);
        my @codonList = ($geneListRef->[$geneIndex{$id}]->seq =~ m/.{3}/g);
        my $alnseq;
        while(my $aa = shift(@aaList)){
            $alnseq .= ($aa eq "-") ? "---" : shift(@codonList);
        }
        $alnseq .= shift(@codonList) if(@codonList);
        $geneListRef->[$geneIndex{$id}]->seq($alnseq); #Keeps the Actual Sequence
        $alnseq =~ tr/Nn/--/; #Treats N as uniformative for distance calculations
        warn "[WARNING] Unused Codons in alignment of $id for $geneName\n" if(@codonList);
        $seqObj->seq($alnseq);
    }
    #unless(close($fh)){
    #    warn "[WARNING] Aligner did not execute correctly for $geneName\n";
    #}
    return $alnObj;
}

#Given a set of aligned sequences outputs a neighbor joining tree
#Input: a Bio::Align::AlignI compliant object, and a genename
#Output: a Bio::Tree Object
sub ConstructTree($$){
    my ($alnObj,$geneName) = @_;
    my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
    my $stats = Bio::Align::DNAStatistics->new;
    my $matObj = $stats->distance(-method => "Uncorrected", -align => $alnObj);
    my $treeObj = $dfactory->make_tree($matObj);
    ReRootTree($treeObj);
    return $treeObj;
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

#Given a set of genes and a tree describing their relationships, extract the representative sequences
# for each cluster
#Input: a reference to an array of GENE objects, a Bio::Tree Obj, and a gene name
#Output: a array of CLUSTER objects
sub ExtractClusters($$$$){   #Usage: my @clusterList = ExtractClusters(\@geneList,$treeObj,$geneName);
    my ($geneListRef,$alnObj,$treeObj,$geneName) = @_;
    my $nextID = 0;
    ##Determine to which cluster each tip belongs
    my @ClusteredLabels = PartitionTipLabels($treeObj->get_root_node,$opts{d});
    my @clusterList;
    ##Construct the Clusters
    foreach my $cluster (@ClusteredLabels){
        my $clustID = sprintf("%s_%0*d",$geneName,$CLUST_ID_PADDING,$nextID++);
        ##Create a New Cluster, with the consensus sequence
        my $consensus = $alnObj->select_noncont_by_name(@{$cluster})->consensus_string();
        $consensus =~ s/\?//g;
        my $clustObj = CLUSTER->new(uid => $clustID, seq => $consensus);
        foreach my $Label (sort @{$cluster}){
            ##Search the seq list for the label and move the GENE obj to the Cluster
            my $i = 0;
            $i++ until($geneListRef->[$i]->uid eq $Label);
            push(@{$clustObj->members},splice(@{$geneListRef},$i,1));
        }
        ##Add Cluster to the List
        push(@clusterList,$clustObj);
    }
    return @clusterList;
}

#Recursive method for partitioning tips of a tree into clusters which have a maximum internal root-tip
#distance
#Input: a Bio::Tree::Node object in a tree, and a max root-tip distance threshold
#Output: an array of array references, where each array ref is tip labels for one cluster
sub PartitionTipLabels($$){
    my ($node,$thresh) = @_;
    if($node->is_Leaf){
        return [$node->id];
    }
    if($node->height <= $thresh){
        my @leaves = $node->get_all_Descendents;
        for(my $i = 0; $i < @leaves; $i++){
            unless($leaves[$i]->is_Leaf){
                splice(@leaves,$i,1);
                $i--;
            }
        }
        return ([map {$_->id} @leaves]);
    }
    return  map {PartitionTipLabels($_,$thresh)} $node->each_Descendent;
}
