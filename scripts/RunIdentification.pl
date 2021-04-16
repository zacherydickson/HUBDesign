#!/usr/bin/perl
package CLUSTER;
use Class::Struct;

struct (CLUSTER => {id => '$', seq =>'$'});

package PSEUDORG;
use Class::Struct;
use Scalar::Util qw(looks_like_number);

struct (PSEUDORG => {id => '$'});
#Has readonly members: len, seq, count
#Has private members: clustList, posList
#   clustList is a list of CLUSTER Objects
#   posList is the position within the buffered pseudo-organism sequence at which each cluster object
#   begins
#Has functions add, map_pos, new

#Initializes a PSEUDORG obj, alternate to new
sub init (){
    my $self = shift;
    if(@_){
        my %Args = @_;
        $self->id($Args{id});
        delete($Args{id});
        warn("[WARNING] too many arguments to PSEUDORG::init\n") if(scalar(keys %Args));
    }
    $self->{'PSEUDORG::len'} = 0;
    $self->{'PSEUDORG::clustList'} = [];
    $self->{'PSEUDORG::posList'} = [];
    return $self;
}

#Adds a cluster object to the Pseudo-organism
#Returns 0 if unsuccessful, and the current number of clusters otherwise
sub add (){
    my ($self,$clustObj) = @_;
    unless(defined $clustObj){
        warn "[WARNING] Not enough arguments to PSEUDORG::add\n";
        return 0;
    }
    unless(ref($clustObj) == "CLUSTER"){
        warn "[WARNING] Cannot PSEUDORG::add non-cluster to PSEUDORG $self->id\n";
        return 0;
    }
    unless(defined $clustObj->id and defined $clustObj->seq){
        warn "[WARNING] Attempt to add at least partially undefined cluster to PSEUDORG $self->id\n";
        return 0;
    }
    push(@{$self->{'PSEUDORG::clustList'}},$clustObj);
    push(@{$self->{'PSEUDORG::posList'}}, $self->{'PSEUDORG::len'} + 1);
    $self->{'PSEUDORG::len'} += length($clustObj->seq) + 1;
    return scalar(@{$self->{'PSEUDORG::clustList'}});
}

#given a 1 indexed position in the buffered pseudo-organism sequence
#returns a two element vector, the id of the cluster that position is in, and the position within that
#cluster
#returns zero if unsuccessful
sub map_pos(){
    my ($self,$pos) = @_;
    unless(defined $pos){
        warn "[WARNING] Not enough arguments to PSEUFORG::map_pos\n";
        return 0;
    }
    unless(looks_like_number($pos) and $pos >= 1 and $pos <= $self->len){
        warn "[WARNING] Cannot PSEUDORG::map_pos which is outside the sequence\n";
        return 0;
    }
    my $index = 0;
    $index++ while($index + 1 < $self->count and $self->{'PSEUDORG::posList'}->[$index + 1] <= $pos);
    $pos = $pos - $self->{'PSEUDORG::posList'}->[$index] + 1;
    #Check if the mapped position happens to the the buffer position between sequences;
    if($pos > length($self->{'PSEUDORG::clustList'}->[$index]->seq)){
        $pos = undef;
    } 
    return ($self->{'PSEUDORG::clustList'}->[$index]->id, $pos);
}

#Returns the length of the buffered pseudo-organism sequence
sub len(){
    my $self = shift;
    return $self->{'PSEUDORG::len'};
}

#Returns the buffered pseudo-organism sequence
sub seq(){
    my $self = shift;
    return join("n", (map { uc($_->seq) } @{$self->{'PSEUDORG::clustList'}}));
}

#Return sthe number of clusters in the psuedo-organism
sub count(){
    my $self = shift;
    return scalar(@{$self->{'PSEUDORG::clustList'}});
}

sub get_ids(){
    my $self = shift;
    return map {$_->id} @{$self->{'PSEUDORG::clustList'}};
}

package BAITREGION;
use Class::Struct;

struct (BAITREGION => {taxon_id => '$', clust_id => '$', pos => '$', seq => '$'});

sub isContig(){
    my ($self,$other) = @_;
    unless(ref($self) eq "BAITREGION" and ref($other) eq "BAITREGION"){
        warn "[WARNING] Cannot check contiguity of non BAITREGION objects\n";
        return 0;
    }

    if($self->taxon_id eq $other->taxon_id){
        if($self->clust_id eq $other->clust_id){
            ($self,$other) = ($other,$self) if($other->pos < $self->pos);
            if($self->pos + length($self->seq) == $other->pos + length($other->seq) - 1){
                return 1;
            }
        }
    }
    return 0;
}

sub toStr(){
    my $self = shift;
    return join("\t",($self->taxon_id,$self->clust_id,$self->pos,$self->seq));
}

package main;
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use Bio::SeqIO;
use File::Temp;

#============================================================
#Declarations

my $_BOND_EXECUTABLE = "SA_BOND";
my $_ID_PADDING = 6;
my $MIN_PERC = 0;
my $MAX_PERC = 100;
my $MIN_PROBE_LENGTH = 1;
my $MIN_TM_RANGE = 1;
my %DEFAULT;
@DEFAULT{qw(d g i l m n p r t k)} = ('/tmp', '30.0:70.0', 75 ,75, '10.0', -1,'50.0', 15, 0, 0);
my %LOG_LEVEL = (ERROR => 0, WARNING => 1, INFO => 2, DEBUG => 3);

sub GetProcessorCount();	 #my $cpu_count = GetProcessorCount();
sub LogParameters();             #LogParameters();
sub LogEvent($$);		 #LogEvent($message,$level)
sub ProcessNumericOption($$$$$$);#my $opt = ProcessNumericOption($val,$default,$min,$max,$bInt,$varName)
sub OpenFileHandle($$;$);	 #my $filehandle = OpenFileHandle($file,$file_type,\&exit_func);
sub ParseAssignments($);         #my %ClusterAssignments = ParseAssignments($assignmentFile);
sub ConstructPseudoGenomes($$);	 #my %PGDict = ConstructPseudogenomes($clustersSeqFile,\%ClustAssign);
sub WritePseudoGenomes($);	 #my $tempFile = WritePseudoGenomes(\%PGDict);
sub RunBOND($);			 #my $tempFile = RunBOND($PGSeqFile);
sub CollapseBaits($$);	         #my @BaitRegions = CollapseBaits($oligFile,\%PGDict);
sub OutputResults($);		 #OutputResults(\@BRList);

#============================================================
#Handling Input

my %opts; #Global Hash with all options
getopts('d:g:i:l:m:n:p:r:t:kvh',\%opts);

if(@ARGV < 2 or exists $opts{h}){
    my $Usage = basename($0). "Assignment.tsv ClustSeq.fna > CandidateBaitRegions.tsv";
    if(exists $opts{h}){
        print STDERR "\n===Description\n".
            "Given a set of representative sequences and information about their heirarchical assignments\n".
            "\tIdentifies candidate probe sequences which are sufficently specific\n". 
            "===Usage\n$Usage\n".
            "===Inputs\n".
            "Assignment.tsv\tA tab separated values with headers of ClusterID, TaxonID, and Penetrance\n".
            "ClustSeq.fna\tA fasta formated file with headers which match ClusterIDs in Assignment.tsv\n".
            "===Outputs\n".
            "A tab separated file with headers of RegionID, TaxonID, ClustID, Position, and Sequence\n".
            "===Options\n".
            "-d dir\tiDirectory for temporary pseudo-genome files (Default $DEFAULT{d})\n".
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
            "===Flags\n".
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
#Handle TmpDir
unless (-d $opts{d}){
    LogEvent("Directory $opts{d} does not exists","ERROR");
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
    LogEvent("GC parameters must be in the format min:max","ERROR");
}
$opts{g} = [split(/:/,$opts{g})];
for(my $i = 0; $i < 2; $i++){
    my $val = $opts{g}->[$i];
    if($val eq ""){
        my $boundary = $MIN_PERC;
        $boundary = $MAX_PERC if($i == 1);
        $opts{g}->[$i] = $boundary;
    } elsif(!looks_like_number($val) or $val < $MIN_PERC or $val > $MAX_PERC){
        LogEvent("GC Parameters must be bewtween 0 and 100, inclusive","ERROR");
    }
}
@opts{qw(g_min g_max)} = @{$opts{g}};
unless(defined $opts{g_min} or defined $opts{g_max}){
    LogEvent("At least one of min or max must be specified for GC parameters","ERROR");
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
            LogEvent("Melting Temperature extrema must be valid celsius temperatures","ERROR");
        }
    }
    @opts{qw(m_min m_max)} = @{$opts{m}};
    unless(defined $opts{m_min} or defined $opts{m_max}){
        LogEvent("At least one of min or max must be specified for T_m range","ERROR");
    }
    if(defined $opts{m_min} and defined $opts{m_max} and $opts{m_min} > $opts{m_max}){
        @opts{qw(m_min m_max)} = @opts{qw(m_max m_min)};
    }
} else { #Interval Size specified
    $opts{m_range} = ProcessNumericOption($opts{m},$DEFAULT{m},1,undef,0,"T_m Interval Size");
    
}
#Ensure Max Oligs is non-zero
if($opts{n} == 0){
    LogEvent("Max Probes must be non-zero","ERROR");
}
#Ensure Threads is within system limits
my $cpu_count = GetProcessorCount();
if($cpu_count == -1){
    LogEvent("Could not determine Sys_max threads: $!\n\t Can't fully validate max_threads\n","WARNING");
    if($opts{t} == 0){
        LogEvent("Sys_max threads unknown: proceeding with 1 thread","WARNING");
        $opts{t} = 1;
    }
} elsif($opts{t} == 0 or $opts{t} > $cpu_count){
    if($opts{t} > $cpu_count){
        LogEvent("Max threads is greater than sys_max: proceeding with $cpu_count threads","WARNING");
    }
    $opts{t} = $cpu_count;
}
#Handle Verbosity
$opts{v} = (exists $opts{v}) ? 2 : 1;

#Handle Arguments
my %FileDict;
@FileDict{qw(Assignment ClustSeq)} = @ARGV[0 .. 1];

while( my ($type,$file) = each %FileDict){
    LogEvent("$type file ($file) does not exists","ERROR") unless (-e $file);
}

#============================================================
#Main Script

LogParameters();
my %ClusterAssignment = ParseAssignments($FileDict{Assignment});
my %PseudoGenomeDict = ConstructPseudoGenomes($FileDict{ClustSeq},\%ClusterAssignment);
$FileDict{PseudoGenome} = WritePseudoGenomes(\%PseudoGenomeDict);
$FileDict{Oligos} = RunBOND($FileDict{PseudoGenome});
my @BaitRegionList = CollapseBaits($FileDict{Oligos},\%PseudoGenomeDict);
OutputResults(\@BaitRegionList);

#============================================================
#Subroutines

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

#Outputs the Parameters to STDERR
sub LogParameters(){
    return if($opts{v} < $LOG_LEVEL{INFO});
    select(STDERR);
    print join("",(('=') x 60)),"\n";
    print "Initialized on ".(localtime)."\n";
    print "Parameters:\n";
    my @paramNames = ("Probe Length","Max Similarity","Max Run","GC Range","T_m Interval","Min Penetrance","Threads","Assignments","Clusters","Temp Directory","Keep Files");
    my %params;
    @params{@paramNames} = ($opts{l},$opts{i},$opts{r},"$opts{g_min} - $opts{g_max}",undef,$opts{p},$opts{t},$FileDict{Assignment},$FileDict{ClustSeq},$opts{d},undef);
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
    foreach my $param (@paramNames){
        printf("%15s : %s\n",$param,$params{$param});
    }
    print join("",(('=') x 60)),"\n";
    select(STDOUT);
}

#Writes a message to STDERR with a timestamp and the level of severity, exits if necessary
sub LogEvent($$){
    my ($message,$level) = @_;
    return if($LOG_LEVEL{$level} > $opts{v});
    my ($sec,$min,$hour) = localtime;
    select(STDERR);
    printf("%02d:%02d:%02d - [%s] %s\n",$hour,$min,$sec,$level,$message);
    exit 1 if($LOG_LEVEL{$level} == 0);
    select(STDOUT);
}

sub ProcessNumericOption($$$$$$){
    my($val,$default,$min,$max,$bInt,$varName) = @_;
    return $default unless(defined $val);
    if(looks_like_number($val)){
        $val = int($val) if($bInt);
        if(!defined $min or $val >= $min){
            if(!defined $max or $val <= $max){
                return $val;
            }
        }
    }
    my $message = sprintf("%s must be a%s between %s and %s",$varName,
       ($bInt ? "n integer" : " value"),(defined $min ? $min : "-∞"),(defined $max ? $max : "∞"));
    LogEvent($message,"ERROR");
}

#Given a file path, opens it for reading and fails with the provided exit function if unsuccessful
# The default exit function is die
#Input	File - Path to the file to open
#	Type - A descriptor of the type of file being opened
# (opt)	exit_Func - a reference to a subroutine to run if opening the file fails
#Output	A file handle if successfule, 0 otherwise.
#Note: Calling process must close the file handle
sub OpenFileHandle($$;$){
    my ($file,$type,$level) = @_;
    $level = "ERROR" unless(defined $level);
    if(open(my $fh, $file) ){
        return $fh;
    } else {
        LogEvent("Could not open $type file ($file): $!",$level);
        return 0;
    }
}

#Reads a cluster assignment file and retains clusters with sufficient penetration
#Input	File - Path to Cluster Assignment File
#Output	A hash with keys of cluster_id and values of taxon_id
sub ParseAssignments($){
    LogEvent("Parsing cluster assignments...","INFO");
    my $fh = OpenFileHandle(shift,"Assignment");
    my @headers = split(/\t/,<$fh>);
    if(@headers < 3){
        LogEvent("Incorrectly formatted Cluster assignment file: Use -h to see correct format","ERROR");
    }
    my %ClustAssign;
    my $clustCount = 0;
    while(my $line = <$fh>){
        chomp($line);
        my ($clust_id,$taxon_id,$penetrance) = split(/\t/,$line);
        unless(looks_like_number($penetrance) and $penetrance >= 0 and $penetrance <= 100){
            LogEvent("Penetrance for $clust_id not in the range [0,100] - Skipping","WARNING");
            next;
        }
        $clustCount++;
        if($penetrance >= $opts{p}){
            if(exists $ClustAssign{$clust_id}){
                LogEvent("Multiple assignments for $clust_id found: ignoring assignment to $taxon_id","WARNING");
                next;
            }
            $ClustAssign{$clust_id} = $taxon_id;
        }
    }
    close($fh);
    my $passCount = scalar(keys %ClustAssign);
    LogEvent(sprintf("Parsed %d assignments with %d (%0.1f%%) with sufficient penetrance",$clustCount,$passCount,$passCount/$clustCount*100),"INFO");
    return %ClustAssign;
}

sub ConstructPseudoGenomes($$){
    LogEvent("Constructing pseudo-genomes...","INFO");
    my $SeqIOObj = Bio::SeqIO->new(-file => shift);
    my $ClustAssignRef = shift;
    my %PGDict;
    my $clustCount = 0;
    my %SeenClusters;
    while( my $seqObj = $SeqIOObj->next_seq){
        my $clust_id = $seqObj->id;
        next unless(exists $ClustAssignRef->{$clust_id});
        if(exists $SeenClusters{$clust_id}){
            LogEvent("Ignoring duplicate sequnce for $clust_id","WARNING");
            next;
        }
        $SeenClusters{$clust_id} = 1;
        $clustCount++;
        my $taxon_id = $ClustAssignRef->{$clust_id};
        my $clustObj = CLUSTER->new(id => $clust_id, seq => $seqObj->seq);
        my $PG = new PSEUDORG;
        $PGDict{$taxon_id} = $PG->init(id => $taxon_id) unless(exists $PGDict{$taxon_id});
        $PGDict{$taxon_id}->add($clustObj);
    }
    my $taxonCount = scalar(keys %PGDict);
    LogEvent(sprintf("Partitioned %d clusters into %d pseudo-genomes",$clustCount,$taxonCount),"INFO");
    return %PGDict;
}

sub WritePseudoGenomes($){
    my $PGDictRef = shift;
    LogEvent("Writing psuedo-genomes to file...","INFO");
    my $tmpfile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDESIGN_PG_XXXX", SUFFIX => ".fna", UNLINK => ($opts{k}+ 1) % 2);
    my $SeqIOObj = Bio::SeqIO->new(-file => ">$tmpfile",-format => 'fasta');
    my $totalLength = 0;
    foreach my $taxon_id (sort keys %{$PGDictRef}){
        my $seqObj = Bio::Seq->new(-id => $taxon_id, -seq => $PGDictRef->{$taxon_id}->seq);
        $totalLength += $PGDictRef->{$taxon_id}->len;
        $SeqIOObj->write_seq($seqObj);
    }
    LogEvent("Pseudo-genomes totalling $totalLength bp written to $tmpfile","INFO");
    return $tmpfile;
}

#Given a file containing Pseudo genome sequences, Calls BOND and returns the file containing the raw
#bond output
sub RunBOND($){
    my $infile = shift;
    LogEvent("Identifying candidate probes...","INFO");
    my $outfile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDESIGN_CP_XXXX", SUFFIX => ".olig", UNLINK => ($opts{k}+ 1) % 2);
    my $TmFlag;
    if (exists $opts{m_range}) {
        $TmFlag = "-rangeTm $opts{m_range}";
    } else {
        $TmFlag = "-minGC $opts{m_min} -maxGC $opts{m_max}";
    }
    my $cmd = "$_BOND_EXECUTABLE $infile $outfile -length $opts{l} -seqSim $opts{i} -maxMatch $opts{r} $TmFlag -minGC $opts{g_min} -maxGC $opts{g_max} -maxOligs $opts{n} 1>&2";
    LogEvent("Running BOND with command: $cmd","INFO");
    my $result = system($cmd);
    LogEvent("BOND failure: $!","ERROR") unless($result == 0);
    LogEvent("Raw BOND output in $outfile","INFO");
    return $outfile;
}

sub CollapseBaits($$){
    my ($oligFile,$PGDictRef) = @_;
    LogEvent("Collapsing contiguous candidate probes...","INFO");
    my $fh = OpenFileHandle($oligFile,"Raw Oligos","WARNING");
    if(!$fh){
        $oligFile->unlink_on_destroy(0);
        LogEvent("Failure to parse Raw Oligos, Saved to $oligFile", "ERROR");
    }
    my $candidate_probe_count = 0;
    my @BRList;
    while(my $line = <$fh>){
        chomp($line);
        my ($tgtStr,$seq,undef,undef,$posStr,undef) = split(/\t/,$line);
        my (undef,$index) = split(/:/,$tgtStr);
        $index--;
        my (undef,$pos) = split(/:/,$posStr);
        $pos++;
        my $taxon_id = (sort keys %{$PGDictRef})[$index];
        my ($clust_id,$clust_pos) = $PGDictRef->{$taxon_id}->map_pos($pos);
        $candidate_probe_count++;
        my $brObj = BAITREGION->new(taxon_id => $taxon_id,clust_id => $clust_id,pos => $clust_pos, seq => $seq);
        push(@BRList,$brObj);
    }
    LogEvent("Read $candidate_probe_count candidate probes","INFO");
    @BRList = sort {$a->taxon_id cmp $b->taxon_id or $a->clust_id cmp $b->clust_id or $a->pos <=> $b->pos} @BRList;
    #Collapse contiguous probes together
    for(my $i = 1; $i < @BRList; $i++){
        if($BRList[$i-1]->isContig($BRList[$i])){
            $BRList[$i-1]->seq($BRList[$i-1]->seq.substr($BRList[$i]->seq,-1));
            splice(@BRList,$i,1);
            $i--;
        }
    }
    LogEvent("Collapsed probes into ".scalar(@BRList)." candidate bait regions","INFO");
    return @BRList;
}

sub OutputResults($){
    my $BRListRef = shift;
    LogEvent("Outputting Bait Regions...","INFO");
    for(my $i = 0; $i < @{$BRListRef}; $i++){
        my $reg_id = sprintf("BR_%0*d",$_ID_PADDING,$i);
        print  "$reg_id\t",$BRListRef->[$i]->toStr,"\n";
    }
    LogEvent("Done","INFO");
}
