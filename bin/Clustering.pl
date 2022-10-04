#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
#use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;
use Getopt::Std;
use Class::Struct;
use File::Temp;
use FindBin;
use File::Spec;
use lib File::Spec->catdir($FindBin::RealBin, '..', 'lib');
use HUBDesign::Util qw(ValidateThreadCount ProcessNumericOption OpenFileHandle LoadConfig);
use Bio::Align::ParallelDNAStatistics;


#============================================================
#Declarations

struct  (GENE => {tid => '$', chr => '$', uid => '$', start => '$', end => '$', strand => '$', seq => '$'});
struct  (CLUSTER => {uid => '$', seq => '$', members => '@',cigar => '%'});

sub LoadGFFSet($@);         #Usage: my %GFFSet = LoadGFFSet($GFFListFile,@GFFFileList);
sub LoadGeneInfo(%);        #Usage: my %GeneInfo = LoadGenes(%GFFSet);
sub HandleSinglets($$$);    #Usage: my $nSinglets = HandleSinglets(\%GeneInfo,$ClustHandle,$SeqOObj);
sub OutputCluster($$$);     #Usage: OutputCluster($ClustObj,$ClustHandle,$SeqOObj);
sub AlignByProt($$);        #Usage: AlignByProt(\@geneList,$geneName);
sub AlignByDNA($$);         #Usage: AlignByDNA(\@geneList,$geneName);
sub ConstructTree($$);      #Usage: my $TreeObj = ConstructTree($alnObj,$geneName);
sub ReRootTree($);          #Usage: ReRootTree($treeObj);
sub GetCIGARHash($);	    #Usage: my %cigar = GetCIGARHash($alnObj);
sub ExtractClusters($$$$);  #Usage: my @clusterList = ExtractClusters(\@geneList,$alnObj,$treeObj,$geneName);
sub PartitionTipLabels($$); #Usage: my @ClusteredLAbels = PartitionTipLabels($nodeObj,$threshold)
sub ParseConfig($);         #Usage: ParseConfig($ConfigFile);

my $CLUST_ID_PADDING = 5;
my $TEMP_DIR = "/tmp";
my $ALIGNER_EXEC = "mafft";
my $ALIGNER = "$ALIGNER_EXEC --auto --localpair --quiet --thread {threads} {input}";
my $ALIGNER_FORMAT = "fasta";
my %DEFAULT = (d => 0.15, e => '.gz,.gff', p => 1, o => 'ClustInfo.csv', t => 0);

#============================================================
#Handling Input


my %opts;
getopts('d:e:l:p:o:t:C:vVh',\%opts);
ParseConfig($opts{C});

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
            "-d [0-1]\tThe Maximum root-tip divergence within clusters (Default: $DEFAULT{d})\n".
            "-e STR\tA comma separated list of file extensions [Default: $DEFAULT{e}]\n".
            "-l PATH\tA tab delimited file with a path to a gff file, and\n".
            "\t\tan optional taxonID on each line\n".
            "\t\tNote: If -l is specified, then no command line GFF's need to be specified\n".
            "-p INT\tTranslation Table to use (Default: $DEFAULT{p})\n".
            "\t\tNote: Setting translation table to 0 will preventing aligning by protein sequences\n".
            "-o PATH\t The file to which cluster info will be written (Default: $DEFAULT{o})\n".
            "-t STR\tMaximum number of threads, 0 is sys_max (Default: $DEFAULT{t})\n".
            "-C PATH\t Path to a HUBDesign Config file (Default: $FindBin::RealBin/../HUBDesign.cfg)\n".
            "===Flags\n".
            "-v\tVerbose Output\n".
            "-h\tDisplay this message an exits\n";
    } elsif(! exists $opts{l}) {
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

foreach my $opt (keys %DEFAULT){
    $opts{$opt} = $DEFAULT{$opt} unless(exists $opts{$opt});
}

#Process Simple Numeric Options
$opts{d} = ProcessNumericOption($opts{d},$DEFAULT{d},0,1,0,"r2t Divergence");
$opts{p} = ProcessNumericOption($opts{p},$DEFAULT{p},0,33,1,"Translation Table");
$opts{t} = ProcessNumericOption($opts{t},$DEFAULT{t},0,undef,1,"Max Threads");

$opts{d} = (exists $opts{d}) ? $opts{d} : 0.15;
if($opts{d} < 0 or $opts{d} > 1){
    die "Clustering threshold must be a floating value between 0 and 1 inclusive\n";
}
$opts{e} = [split(",",$opts{e})];

if(exists $opts{l}){
    die "Could not locate file: $opts{l}\n" unless(-e $opts{l});
} else{
    $opts{l} = undef;
}

#Handle Num Threads
$opts{t} = ValidateThreadCount($opts{t});
$ALIGNER =~ s/\{threads\}/$opts{t}/;

#Handle Verbosity
$opts{v} = 1 if (exists $opts{V});

#============================================================
#Main Script

my $SeqOObj = Bio::SeqIO->newFh(-format => 'fasta');
open (my $clustHandle, ">$opts{o}") or die "Could not write to $opts{o}: $!\n";
print $clustHandle join("\t",qw(ClusterID TaxonID SeqID GeneID Start End Strand CIGAR)),"\n";


print STDERR "Loading Gene Info...\n" if(exists $opts{v});
my %GeneInfo = LoadGeneInfo(LoadGFFSet($opts{l},@ARGV));
printf STDERR "Gene families Loaded:\t%d\n",scalar(keys %GeneInfo) if(exists $opts{v});
my $nSinglets = HandleSinglets(\%GeneInfo,$clustHandle,$SeqOObj);
printf STDERR "Single member gene families:\t%d\n",$nSinglets if(exists $opts{v});

#Global Declaration of BioPerl Objects for Clustering
my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'UPGMA');
my $statObj;
if($opts{t} == 1){
    $statObj = Bio::Align::DNAStatistics->new();
} else {
    $statObj = Bio::Align::ParallelDNAStatistics->new(-threads => $opts{t});
}

foreach my $geneName (keys %GeneInfo){
    printf STDERR "Processing gene family %s...\n",$geneName if(exists $opts{v});
    my $seqListRef = $GeneInfo{$geneName};
    print STDERR "\tAligning...\n" if(exists $opts{v});
    my $alnObj;
    if($opts{p} != 0){
        $alnObj = AlignByProt($seqListRef,$geneName);
    } else {
        $alnObj = AlignByDNA($seqListRef,$geneName);
    }
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
    print STDERR "|>" if(exists $opts{v});
    my %GeneInfo;
    ##Process Each GFF File
    my $total = scalar keys %gffSet;
    my $count = 0;
    my $next = 1;
    while( my ($tid,$filename) = each %gffSet){
        $count++;
        until($next > int(($count / $total)*100)){
            $next++;
            print STDERR "\b=>" if (exists $opts{v});
        }
        #Handle potentially gzipped files
        my $opencmd = $filename;
        if ($filename =~ m/\.gz$/){
            $opencmd = "gzip -dc $filename |"
        }
        if(open(my $fh, $opencmd)){
            ##Parse the portion of the GFF with gene info
            my %chrIndex; #For rapidly filling in sequence info
            my %chrGNameSetDict; #^
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
                $chrGNameSetDict{$seqID} = {} unless(exists $chrGNameSetDict{$seqID});
                $chrGNameSetDict{$seqID}->{$geneName} = 1;
            }
            ##Extract the Sequence Info
            my $SeqI = new Bio::SeqIO(-fh => $fh, -format => 'fasta');
            while(my $seqObj = $SeqI->next_seq){
                my $id = $seqObj->id;
                next unless(exists $chrGNameSetDict{$id});
                foreach my $geneName (keys %{$chrGNameSetDict{$id}}){
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
    print STDERR "\b=|\n" if (exists $opts{v});
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
            my $clustObj = CLUSTER->new(uid => $clustID, seq => $geneListRef->[0]->seq, members => $geneListRef, cigar => {$geneListRef->[0]->uid => length($geneListRef->[0]->seq)."="});
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
        print join("\t",($clustObj->uid,$geneObj->tid,$geneObj->chr,$geneObj->uid,$geneObj->start,$geneObj->end,$geneObj->strand,$clustObj->cigar->{$geneObj->uid})),"\n";
    }
    print $SeqOObj Bio::Seq->new(-seq => $clustObj->seq, -id => $clustObj->uid);
    select($oldFH);
}


##TODO: Combine Alignment code into a single align function
##TODO: Split De-translation out into its own function

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
    my $cmd = $ALIGNER;
    my $filename = $tmpFile->filename;
    $cmd =~ s/\{input\}/$filename/;
    open( my $fh, "-|", "$cmd") or warn "Alignment Error for $geneName: $!";
    unless($fh){
        die "Alignment error for $geneName: $!";
    }
    my $AlnI = Bio::AlignIO->new(-fh => $fh, -format => $ALIGNER_FORMAT);
    my $alnObj = $AlnI->next_aln();
    unless($alnObj){
        die "Alignment error for $geneName: $!";
    }
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

#Aligns a set of sequences by their dna sequence
#Input: a reference to a list of GENE Objects
#Output: a Bio::Align::AlignI compliant object
sub AlignByDNA($$){
    my ($geneListRef,$geneName) = @_;
    my @Alignment;
    ## Output Translated Sequences to a temporary File
    my $tmpFile = File::Temp->new(DIR => $TEMP_DIR, TEMPLATE => "HUBCLUSTXXXX", SUFFIX =>".fna",UNLINK=> ! exists $opts{V});
    print STDERR "\t\tOutputing sequences to $tmpFile...\n" if(exists $opts{V});
    my $SeqO = Bio::SeqIO->newFh(-format => "fasta", -fh => $tmpFile);
    my %geneIndex;
    my $nextIndex=0;
    foreach my $geneObj (@{$geneListRef}){
        print $SeqO Bio::Seq->new(-seq => $geneObj->seq, -id => $geneObj->uid);
        $geneIndex{$geneObj->uid} = $nextIndex++;
    }
    ## Call the external aligner, and get its output
    my $cmd = $ALIGNER;
    my $filename = $tmpFile->filename;
    $cmd =~ s/\{input\}/$filename/;
    $cmd =~ s/--amino/--nuc/;
    print STDERR "\t\tAligning sequences with cmd\n\t$cmd\n" if(exists $opts{V});
    open( my $fh, "-|", "$cmd") or warn "Alignment Error for $geneName: $!";
    unless($fh){
        die "Alignment error for $geneName: $!";
    }
    my $AlnI = Bio::AlignIO->new(-fh => $fh, -format => $ALIGNER_FORMAT);
    my $alnObj = $AlnI->next_aln();
    unless($alnObj){
        die "Alignment error for $geneName: $!";
    }
    return $alnObj;
}

#Given a set of aligned sequences outputs a neighbor joining tree
#Input: a Bio::Align::AlignI compliant object, and a genename
#Output: a Bio::Tree Object
sub ConstructTree($$){
    my ($alnObj,$geneName) = @_;
    my $matObj = $statObj->distance(-method => "Uncorrected", -align => $alnObj);
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

#Given an alignment returns the CIGAR strings for each sequence with the consensus as a reference
#Input: a Bio::SimpleAlign object
#Output: a hash with keys of seqid and values of SAM style CIGAR strings
sub GetCIGARHash($){
    my $self = shift;
    die "GetCIGARHash requires a Bio::SimpleAlign object" unless(ref($self) eq "Bio::SimpleAlign");
    my @consChars = split(//,$self->consensus_string());
    my %cigar;
    foreach my $seqObj ($self->each_seq){
        my $str;
        my @qChars = split(//,$seqObj->seq);
        my $run = 0;
        my $lastChar = 0;
        my $gapChar = $self->gap_char;
        for(my $i = 0; $i < $self->length; $i++){
            my $curChar;
            if($qChars[$i] eq $consChars[$i]){ #Both the same
                if($qChars[$i] eq $gapChar){ #Both Gaps, skip
                    $curChar = undef;
                } else { #Match
                    $curChar = '=';
                }
            } else { #Both different
                if($qChars[$i] eq $gapChar){ #Deletion
                    $curChar = 'D';
                } elsif($consChars[$i] eq $gapChar){ #Insertion
                    $curChar = 'I';
                } else { #Mismatch
                    $curChar = 'X';
                }
            }
            if(defined $curChar){
                if($curChar eq $lastChar){
                    $run++;
                } else {
                    $str .= $run.$lastChar if($run);
                    $run = 1;
                    $lastChar = $curChar;
                }
            }
        }
        $str .= $run.$lastChar;
        $cigar{$seqObj->id} = $str;
    }
    return %cigar;
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
        my $subalnObj = $alnObj->select_noncont_by_name(@{$cluster});
        my $consensus = $subalnObj->consensus_string();
        $consensus =~ s/\?//g;
        my $clustObj = CLUSTER->new(uid => $clustID, seq => $consensus, cigar => {GetCIGARHash($subalnObj)});
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


sub ParseConfig($){
    my $file = shift;
    $file = "$FindBin::RealBin/../HUBDesign.cfg" unless defined $file;
    unless(-e $file){
        warn "Could not find Config file ($file): Revert to Hardcoded defaults\n";
        return;
    }
    my %Config = LoadConfig($file,"WARNING");
    $CLUST_ID_PADDING = $Config{'ID-padding'} if(exists $Config{'ID-padding'});
    $TEMP_DIR = $Config{'cluster-dir'} if(exists $Config{'cluster-dir'});
    $ALIGNER_EXEC = $Config{'Aligner-exec'} if(exists $Config{'Aligner-exec'});
    $ALIGNER = "$ALIGNER_EXEC $Config{'Aligner-cmd'}" if(exists $Config{'Aligner-cmd'});
    $ALIGNER_FORMAT = $Config{'Aligner-format'} if(exists $Config{'Aligner-format'});
    $DEFAULT{d} = $Config{'r2t-divergence'} if(exists $Config{'r2t-divergence'});
    $DEFAULT{p} = $Config{'translation-table'} if(exists $Config{'translation-table'});
    $DEFAULT{t} = $Config{'threads'} if(exists $Config{'threads'});
}

