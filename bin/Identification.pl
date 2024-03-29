#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use Bio::SeqIO;
use File::Basename;
use File::Temp;
use File::Spec;
use FindBin;
use Class::Struct;
use lib File::Spec->catdir($FindBin::RealBin,'..','lib');
use HUBDesign::Util qw(ValidateThreadCount ProcessNumericOption OpenFileHandle LoadConfig);
use HUBDesign::Logger;
use HUBDesign::BaitRegion;
use HUBDesign::Pseudorg;


#============================================================
#Declarations

struct (CLUSTER => {uid => '$', seq =>'$'});

sub GetProcessorCount();	 #my $cpu_count = GetProcessorCount();
sub LogParameters();             #LogParameters();
sub ParseAssignments($);         #my %ClusterAssignments = ParseAssignments($assignmentFile);
sub ConstructPseudoGenomes($$);	 #my %PGDict = ConstructPseudogenomes($clustersSeqFile,\%ClustAssign);
sub WritePseudoGenomes($);	 #my @tempFiles = WritePseudoGenomes(\%PGDict);
sub RunBOND($;$);		 #my $tempFile = RunBOND($PGSeqFile,$LPSeqFile);
sub CollapseBaits($$);	         #my @BaitRegions = CollapseBaits($oligFile,\%PGDict);
sub OutputResults($);		 #OutputResults(\@BRList);
sub ParseConfig($); 	     #Usage: ParseConfig($ConfigFile);

my $_LOWPEN_ID = "_LOWPEN_";
my $_BOND_EXECUTABLE = "SA_BOND";
my $_ID_PADDING = 6;
my $MIN_PERC = 0;
my $MAX_PERC = 100;
my $MIN_PROBE_LENGTH = 1;
my $MIN_TM_RANGE = 1;
my %DEFAULT;
@DEFAULT{qw(d g i l m n p r t P k)} = ('/tmp', '30.0:70.0', 75 ,75, '10.0', -1,'50.0', 15, 0, 0, 0);


#============================================================
#Handling Input

my %opts; #Global Hash with all options
getopts('d:g:i:l:m:n:p:r:t:C:kvh',\%opts);
ParseConfig($opts{C});

if(@ARGV < 2 or exists $opts{h}){
    my $Usage = basename($0). " Assignment.tsv ClustSeq.fna > CandidateBaitRegions.tab";
    if(exists $opts{h}){
        print STDERR "\n===Description\n".
            "Given a set of representative sequences and information about their heirarchical assignments\n".
            "\tIdentifies candidate probe sequences which are sufficently specific\n". 
            "===Usage\n$Usage\n".
            "===Inputs\n".
            "Assignment.tsv\tA tab separated values with headers of ClusterID, TaxonID, and Penetrance\n".
            "ClustSeq.fna\tA fasta formated file with headers which match ClusterIDs in Assignment.tsv\n".
            "===Outputs\n".
            "A tab delimited file with columns of RegionID, TaxonID, ClustID, Position, and Sequence\n".
            "===Options\n".
            "-d dir\tDirectory for temporary pseudo-genome files (Default $DEFAULT{d})\n".
            "-g str\tA string specifying the GC parameters (Default $DEFAULT{g})\n".
            "\t\tFormat min:max; min and max both between 0 and 100 and either may be optional \n".
            "-i [0-100]\tThe maximum percent identity for non-co-targeting probes (Default: $DEFAULT{i})\n".
            "-l [1-∞)\tThe length of probe sequences (Default: $DEFAULT{l})\n".
            "-m str\tA string specifying the melting temperature parameters (Default: $DEFAULT{m})\n".
            "\t\tSpecified either as a single value for the interval size or min:max\n".
            "\t\t\tFor the former, the  T_m of all probes will fall within an interval of that size\n".
            "\t\tmin and max may be any valid celsius tempurature, and either is optional\n".
            "-n {-1,[1,∞)}\tThe maximum number of candidate probes per taxon, -1 is no limit (Default: -1)\n".
            "-p [0-100]\tThe minimum penetrance for inclusion (Default: $DEFAULT{p})\n".
            "-r [1-∞)\tThe maximum number of consectutive identities for non-co-targeting probes (Default: $DEFAULT{r})\n".
            "-t [0-∞)\tThe number of threads to use; 0 is system max (Default: $DEFAULT{t})\n".
            "-C PATH\t Path to a HUBDesign Config file (Default: $FindBin::RealBin/../HUBDesign.cfg)\n".
            "===Flags\n".
            "-P\tDisable blacklisting of sequences too similar to low penetrance sequences\n".
            "-k\tKeep intermediary files\n".
            "-v\tVerbose Output\n".
            "-h\tDisplay this message and exit\n";
    } else {
        print STDERR "Usage: $Usage\n\tUse -h for more info\n";
    }
    exit 1;
}

foreach (keys %DEFAULT){
    $opts{$_} = exists $opts{$_} ? $opts{$_} : $DEFAULT{$_};
}

#Init Logger Verbosity
$opts{v} = (exists $opts{v}) ? "INFO" : "WARNING";
my $Logger = HUBDesign::Logger->new(level => $opts{v});

#Handle TmpDir
unless (-d $opts{d}){
    $Logger->Log("Directory $opts{d} does not exists","ERROR");
}
#Handle Simple Numeric Options 
$opts{i} = ProcessNumericOption($opts{i},$DEFAULT{i},$MIN_PERC,$MAX_PERC,0,"Max Similarity");
$opts{l} = ProcessNumericOption($opts{l},$DEFAULT{l},1,undef,1,"Probe Length");
$opts{n} = ProcessNumericOption($opts{n},$DEFAULT{n},-1,undef,1,"Max Probes");
$opts{p} = ProcessNumericOption($opts{p},$DEFAULT{p},$MIN_PERC,$MAX_PERC,0,"Max Penetrance");
$opts{r} = ProcessNumericOption($opts{r},$DEFAULT{r},1,undef,0,"Max Identiy Run");
$opts{t} = ProcessNumericOption($opts{t},$DEFAULT{t},0,undef,0,"Max Threads");
#Handle GC Param
unless($opts{g} =~ m/^[^:]*:[^:]*$/){
    $Logger->Log("GC parameters must be in the format min:max","ERROR");
}
$opts{g} = [split(/:/,$opts{g})];
for(my $i = 0; $i < 2; $i++){
    my $val = $opts{g}->[$i];
    if($val eq ""){
        my $boundary = $MIN_PERC;
        $boundary = $MAX_PERC if($i == 1);
        $opts{g}->[$i] = $boundary;
    } elsif(!looks_like_number($val) or $val < $MIN_PERC or $val > $MAX_PERC){
        $Logger->Log("GC Parameters must be bewtween 0 and 100, inclusive","ERROR");
    }
}
@opts{qw(g_min g_max)} = @{$opts{g}};
unless(defined $opts{g_min} or defined $opts{g_max}){
    $Logger->Log("At least one of min or max must be specified for GC parameters","ERROR");
}
if(defined $opts{g_min} and defined $opts{g_max} and $opts{g_min} > $opts{g_max}){
    @opts{qw(g_min g_max)} = @opts{qw(g_max g_min)};
}
#Handle Melting Temp
if($opts{m} =~ m/^[^:]*:[^:]*$/){ # Interval limits specified
    $opts{m} = [split(/:/,$opts{m})];
    for(my $i = 0; $i < 2; $i++){
        my $val = $opts{m}->[$i];
        if(!defined $val or $val eq ""){
            $opts{m}->[$i] = undef;
        } elsif(!looks_like_number($val)){
            $Logger->Log("Melting Temperature extrema must be valid celsius temperatures","ERROR");
        }
    }
    @opts{qw(m_min m_max)} = @{$opts{m}};
    unless(defined $opts{m_min} or defined $opts{m_max}){
        $Logger->Log("At least one of min or max must be specified for T_m range","ERROR");
    }
    if(defined $opts{m_min} and defined $opts{m_max} and $opts{m_min} > $opts{m_max}){
        @opts{qw(m_min m_max)} = @opts{qw(m_max m_min)};
    }
} else { #Interval Size specified
    $opts{m_range} = ProcessNumericOption($opts{m},$DEFAULT{m},1,undef,0,"T_m Interval Size");
    
}
#Ensure Max Oligs is non-zero
if($opts{n} == 0){
    $Logger->Log("Max Probes must be non-zero","ERROR");
}
#Ensure Threads is within system limits
$opts{t} = ValidateThreadCount($opts{t});

#Handle Arguments
my %FileDict;
@FileDict{qw(Assignment ClustSeq)} = @ARGV[0 .. 1];

while( my ($type,$file) = each %FileDict){
    $Logger->Log("$type file ($file) does not exist","ERROR") unless (-e $file);
}

#============================================================
#Main Script

LogParameters();
my %ClusterAssignment = ParseAssignments($FileDict{Assignment});
my %PseudoGenomeDict = ConstructPseudoGenomes($FileDict{ClustSeq},\%ClusterAssignment);
@FileDict{qw(PseudoGenome LowPen)} = WritePseudoGenomes(\%PseudoGenomeDict);
$FileDict{Oligos} = RunBOND($FileDict{PseudoGenome},$FileDict{LowPen});
my @BaitRegionList = CollapseBaits($FileDict{Oligos},\%PseudoGenomeDict);
OutputResults(\@BaitRegionList);

#============================================================
#Subroutines

#Outputs the Parameters to STDERR
sub LogParameters(){
     my @paramNames = ("Probe Length","Max Similarity","Max Run","GC Range","T_m Interval","Max Candidates","Min Penetrance","Threads","Assignments","Clusters","Temp Directory","Keep Files","Blacklist LowPen");
    my %params;
    @params{@paramNames} = ($opts{l},$opts{i},$opts{r},"$opts{g_min} - $opts{g_max}",undef,$opts{n}, $opts{p},$opts{t},$FileDict{Assignment},$FileDict{ClustSeq},$opts{d},undef);
    if(exists $opts{m_range}){
        $params{"T_m Interval"} = $opts{m_range};
    } else {
        splice(@paramNames,4,1,"T_m Range");
        delete($params{"T_m Interval"});
        my $min = defined($opts{m_min}) ? $opts{m_min} : "-∞";
        my $max = defined($opts{m_max}) ? $opts{m_max} : '+∞';
        $params{"T_m Range"} = "$min - $max";
    }
    $params{"Keep Files"} = $opts{k} ? "TRUE" : "FALSE";
    $params{"Blacklist LowPen"} = $opts{P} ? "FALSE" : "TRUE";
    $Logger->LogParameters(map {($_, $params{$_})} (sort keys %params));
}

# Reads a cluster assignment file and retains clusters with sufficient penetration
# Inputs -File - Path to Cluster Assignment File
# Output	- A hash with keys of cluster_id and values of taxon_id
sub ParseAssignments($){
    $Logger->Log("Parsing cluster assignments...","INFO");
    my $fh = OpenFileHandle(shift,"Assignment","ERROR");
    my @headers = split(/\t/,<$fh>);
    if(@headers < 3){
        $Logger->Log("Incorrectly formatted Cluster assignment file: Use -h to see correct format","ERROR");
    }
    my %ClustAssign;
    my $clustCount = 0;
    my $failCount = 0;
    while(my $line = <$fh>){
        chomp($line);
        my ($clust_id,$taxon_id,$penetrance) = split(/\t/,$line);
        unless(looks_like_number($penetrance) and $penetrance >= 0 and $penetrance <= 100){
            $Logger->Log("Penetrance for $clust_id not in the range [0,100] - Skipping","WARNING");
            next;
        }
        $clustCount++;
        if($penetrance >= $opts{p}){
            if(exists $ClustAssign{$clust_id}){
                $Logger->Log("Multiple assignments for $clust_id found: ignoring assignment to $taxon_id","WARNING");
                next;
            }
            $ClustAssign{$clust_id} = $taxon_id;
        } elsif (!$opts{P}) {
            $ClustAssign{$clust_id} = $_LOWPEN_ID;
            $failCount++;
        }
    }
    close($fh);
    my $passCount = scalar(keys %ClustAssign) - $failCount;
    $Logger->Log(sprintf("Parsed %d assignments with %d (%0.1f%%) with sufficient penetrance",$clustCount,$passCount,$passCount/$clustCount*100),"INFO");
    return %ClustAssign;
}

# Given a file of cluster sequence and a set of assignments, filters and partitions those clusters into
# pseudo-genomes
# Inputs - A path to a fasta formated file containing cluster sequences
#        - a reference to a hash with keys of cluster ids and values of CLUSTER objects
# Output - a hash with keys of taxon_id and values of HUBDesign::Pseudorg objects
# Note assignment keys and cluster ids are expected to match
sub ConstructPseudoGenomes($$){
    $Logger->Log("Constructing pseudo-genomes...","INFO");
    my $SeqIOObj = Bio::SeqIO->new(-file => shift);
    my $ClustAssignRef = shift;
    my %PGDict;
    my $clustCount = 0;
    my %SeenClusters;
    while( my $seqObj = $SeqIOObj->next_seq){
        my $clust_id = $seqObj->id;
        next unless(exists $ClustAssignRef->{$clust_id});
        if(exists $SeenClusters{$clust_id}){
            $Logger->Log("Ignoring duplicate sequnce for $clust_id","WARNING");
            next;
        }
        $clustCount++;
        $SeenClusters{$clust_id} = 1;
        my $taxon_id = $ClustAssignRef->{$clust_id};
        my $clustObj = CLUSTER->new(uid => $clust_id, seq => $seqObj->seq);
        $PGDict{$taxon_id} = HUBDesign::Pseudorg->new(id => $taxon_id) unless(exists $PGDict{$taxon_id});
        $PGDict{$taxon_id}->add($clustObj);
    }
    my @taxon_idList = keys %PGDict;
    foreach my $taxon_id (@taxon_idList){
        delete $PGDict{$taxon_id} unless($PGDict{$taxon_id}->len >= ($opts{l} +25));
    }
    my $taxonCount = scalar(keys %PGDict);
    $taxonCount-- if(exists $PGDict{$_LOWPEN_ID});
    $Logger->Log(sprintf("Partitioned %d clusters into %d pseudo-genomes; %d pg's were too short",$clustCount,$taxonCount,scalar(@taxon_idList)-$taxonCount),"INFO");
    return %PGDict;
}

#Given a set of pseudo genomes, writes them to a file in fasta format
#Inputs - a reference to a hash of HUBDesign::Pseudorg objects
#Outputs - a File::Temp object to which the sequences were wrtitten
sub WritePseudoGenomes($){
    my $PGDictRef = shift;
    $Logger->Log("Writing psuedo-genomes to file...","INFO");
    my $tmpfile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDESIGN_PG_XXXX", SUFFIX => ".fna", UNLINK => ($opts{k}+ 1) % 2);
    my $SeqIOObj = Bio::SeqIO->new(-file => ">$tmpfile",-format => 'fasta');
    my $lowpenfile = undef;
    my $lowpenIOObj = undef;
    if(!$opts{P} and exists $PGDictRef->{$_LOWPEN_ID}){
        $lowpenfile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDESIGN_LP_XXXX", SUFFIX => ".fna", UNLINK => ($opts{k}+ 1) % 2);
        $lowpenIOObj = Bio::SeqIO->new(-file => ">$lowpenfile",-format => 'fasta');
    }
    my $totalLength = 0;
    foreach my $taxon_id (sort keys %{$PGDictRef}){
        my $seqObj = Bio::Seq->new(-id => $taxon_id, -seq => $PGDictRef->{$taxon_id}->seq);
        if($taxon_id eq $_LOWPEN_ID){
            $lowpenIOObj->write_seq($seqObj);
            $Logger->Log("Low Penetrance sequences totalling @{[$PGDictRef->{$taxon_id}->len]} bp written to $lowpenfile","INFO");
            next;
        }
        $totalLength += $PGDictRef->{$taxon_id}->len;
        $SeqIOObj->write_seq($seqObj);
    }
    $Logger->Log("Pseudo-genomes totalling $totalLength bp written to $tmpfile","INFO");
    return ($tmpfile,$lowpenfile);
}

#Given a file containing Pseudo genome sequences, Calls BOND and returns the file containing the raw
#bond output
#Inputs - A file containing fasta formtatted pseudor-genomes
#Output - A File::Temp object containing the raw BOND output;
sub RunBOND($;$){
    my $infile = shift;
    my $blfile = shift;
    $Logger->Log("Identifying candidate probes...","INFO");
    my $outfile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDESIGN_CP_XXXX", SUFFIX => ".olig", UNLINK => ($opts{k}+ 1) % 2);
    my $TmFlag;
    if (exists $opts{m_range}) {
        $TmFlag = "-rangeTm $opts{m_range}";
    } else {
        $TmFlag = "-minTm $opts{m_min}" if(defined $opts{m_min});
        $TmFlag = join(" ",($TmFlag,"-maxTm $opts{m_max}")) if(defined $opts{m_max});
    }
    my $blFlag = (defined $blfile) ? "-blacklist $blfile" : "";
    #Note: Changes to TilingLen mode affect which PGs are included (ConstructPG)
    my $cmd = "$_BOND_EXECUTABLE $infile $outfile -length ".($opts{l}+25)." -tilingLen $opts{l} -seqSim $opts{i} -maxMatch $opts{r} $TmFlag -minGC $opts{g_min} -maxGC $opts{g_max} -maxOligs $opts{n} $blFlag -maxThread $opts{t} 1>&2";
    $Logger->Log("Running BOND with command: $cmd","INFO");
    $cmd = "stdbuf -i0 -o0 -e0 $cmd" unless(system("stdbuf -e0 sleep 1"));
    my $result = system($cmd);
    $Logger->Log("BOND failure: $result","ERROR") unless($result == 0);
    $Logger->Log("Raw BOND output in $outfile","INFO");
    return $outfile;
}

#Given raw baond output, collpases contiguous probes into bait regions
#Inputs - A path/File:Temp object containing raw BOND output
#       - A reference to a hash of HUBDesign::Pseudorg objects
#Output - A List of HUBDesign::BaitRegion objects
sub CollapseBaits($$){
    my ($oligFile,$PGDictRef) = @_;
    $Logger->Log("Collapsing contiguous candidate probes...","INFO");
    my $fh = OpenFileHandle($oligFile,"Raw Oligos","WARNING");
    if(!$fh){
        $oligFile->unlink_on_destroy(0);
        $Logger->Log("Failure to parse Raw Oligos, Saved to $oligFile", "ERROR");
    }
    my $candidate_probe_count = 0;
    my @BRList;
    my @taxon_idList = (sort keys %{$PGDictRef});
    while(my $line = <$fh>){
        chomp($line);
        my ($tgtStr,$seq,undef,undef,$posStr,undef) = split(/\t/,$line);
        my (undef,$index) = split(/:/,$tgtStr);
        $index--;
        my (undef,$pos) = split(/:/,$posStr);
        $pos++;
        my $taxon_id = $taxon_idList[$index];
        next if($taxon_id eq $_LOWPEN_ID);
        my ($clust_id,$clust_pos) = $PGDictRef->{$taxon_id}->map_pos($pos);
        $candidate_probe_count++;
        my $brObj = HUBDesign::BaitRegion->new(taxon_id => $taxon_id,clust_id => $clust_id,pos => $clust_pos, seq => $seq);
        push(@BRList,$brObj);
    }
    $Logger->Log("Read $candidate_probe_count candidate probes","INFO");
    @BRList = sort {$a->taxon_id cmp $b->taxon_id or $a->clust_id cmp $b->clust_id or $a->pos <=> $b->pos} @BRList;
    #Collapse contiguous probes together
    for(my $i = 1; $i < @BRList; $i++){
        if($BRList[$i-1]->isContig($BRList[$i])){
            $BRList[$i-1]->seq($BRList[$i-1]->seq.substr($BRList[$i]->seq,-1));
            splice(@BRList,$i,1);
            $i--;
        }
    }
    $Logger->Log("Collapsed probes into ".scalar(@BRList)." candidate bait regions","INFO");
    return @BRList;
}

#Given a set of bait regions outputs them to stdout
#Inputs - A reference to an array of HUBDesign::BaitRegion Objects
#Outputs - None, outputs to STDOUT
sub OutputResults($){
    my $BRListRef = shift;
    $Logger->Log("Outputting Bait Regions...","INFO");
    for(my $i = 0; $i < @{$BRListRef}; $i++){
        my $reg_id = sprintf("BR_%0*d",$_ID_PADDING,$i);
        print  "$reg_id\t",$BRListRef->[$i]->toStr,"\n";
    }
    $Logger->Log("Done","INFO");
}

sub ParseConfig($){
    my $file = shift;
    $file = "$FindBin::RealBin/../HUBDesign.cfg" unless defined $file;
    unless(-e $file){
        warn "Could not find Config file ($file): Revert to Hardcoded defaults\n";
        return;
    }
    my %Config = LoadConfig($file,"WARNING");
    $_BOND_EXECUTABLE = $Config{'BOND-exec'} if(exists $Config{'BOND-exec'});
    $_BOND_EXECUTABLE =~ s/\{bin\}/$FindBin::RealBin/;
    $_ID_PADDING = $Config{'ID-padding'} if(exists $Config{'ID-padding'});
    $DEFAULT{t} = $Config{'threads'} if(exists $Config{'threads'});
    $DEFAULT{g} = join(":",((exists $Config{'gc-min'}) ? $Config{'gc-min'} : "", (exists $Config{'gc-max'}) ? $Config{'gc-max'} : ""));
    $DEFAULT{m} = join(":",((exists $Config{'tm-min'}) ? $Config{'tm-min'} : "", (exists $Config{'tm-max'}) ? $Config{'tm-max'} : ""));
    $DEFAULT{m} =~ s/None//g;
    $DEFAULT{m} = $Config{'tm-range'} if($DEFAULT{m} eq ":" and exists $Config{'tm-range'});
    $DEFAULT{i} = $Config{'similarity'} if(exists $Config{'similarity'});
    $DEFAULT{l} = $Config{'length'} if(exists $Config{'length'});
    $DEFAULT{n} = $Config{'candidate-pool'} if(exists $Config{'candidate-pool'});
    $DEFAULT{n} = -1 if($DEFAULT{n} eq 0);
    $DEFAULT{p} = $Config{'penetrance'} if(exists $Config{'penetrance'});
    $DEFAULT{r} = $Config{'contiguity'} if(exists $Config{'match-max'});
}

