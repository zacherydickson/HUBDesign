#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Std;
use File::Basename;
use File::Spec;
use File::Temp;
use lib File::Spec->catdir(
            File::Basename::dirname(File::Spec->rel2abs($0)),
            '..','lib');
use HUBDesign::Util qw(ValidateThreadCount ProcessNumericOption OpenFileHandle);
use HUBDesign::Logger;
use HUBDesign::BaitRegion;

#============================================================
#Declarations

my $_BLAST_EXECUTABLE = "blastn";
my $_LCR_EXECUTABLE = "/home/zac/applications/sdust-master/sdust_";
my $_ID_PADDING = 6;
my $MIN_PERC = 0;
my $MAX_PERC = 100;
my $MIN_PROBE_LENGTH = 1;
my %DEFAULT;
@DEFAULT{qw(d e f i l m t x k)} = ('/tmp', 10 ,'-w 64 -t 20', 75, 75, 20, 0, 'none', 0);

sub LoadBaitRegions($);		#my %BRDict = LoadBaitRegions($brfile)
sub OutputFasta($$);		#my $tmpFile = OutputFasta(\$BRDict,$label);
sub RunLCRFilter($);		#my %ExIVDict = RunLCRFilter($fastafile);
sub ExcludeIntervals($$);	#$ExcludeIntervals(\%BRDict,\%ExIVDict);
sub RunBlastFilter($$);		#my %ExIVDict = RunBlastFilter($fastaFile,$blastdb);
sub OutputResults($);		#OutputResults([values %BRDict]);

#============================================================
#Handling Input

my %opts;
getopts('d:e:f:i:l:m:t:x:kvh',\%opts);

if(@ARGV < 1 or exists $opts{h}){
    my $Usage = basename($0)." [-options] Candidates.tab [BlastDB_1 ...] > FilteredCandidates.tab";
    if(exists $opts{h}){
        print STDERR "\n===Description\n".
            "Given Candidate Bait Regions and filtering parameters, outputs the filtered bait regions\n".
            "===Usage\n$Usage\n".
            "===Required Inputs\n".
            "Candidates.tab\tA tab delimited files with columns of RegionID,TaxonID,ClusterID,Position, and Sequence\n".
            "===Optional Inputs\n".
            "BlastDB_N\t A Path to a blast database to filter against\n".
            "===Outputs\n".
            "A tab delimited files with columns of RegionID,TaxonID,ClusterID,Position, and Sequence\n".
            "===Options\n".
            "-d dir\tiDirectory for temporary fasta files (Default: $DEFAULT{d})\n".
            "-e float\tThe maximum allowable e-value for hits (Default: $DEFAULT{e})\n".
            "-f str\tParameters for the lcr filtering program, use 'no' to disable (Default: $DEFAULT{f})\n".
            "-i [0-100]\tThe maximum allowable percent identity for an hsp (Default: $DEFAULT{i})\n".
            "-l [1-∞)\tThe length of probe sequences (Default: $DEFAULT{l})\n".
            "-m [0-∞)\tThe maximum allowable hsp length (Default: $DEFAULT{m})\n".
            "-t [0-∞)\tThe number of threads to use; 0 is system max (Default: $DEFAULT{t})\n".
            "-x str\tAdditional arguments to the blast cmd (Default: $DEFAULT{x})\n".
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

#Handle Verbosity
$opts{v} = (exists $opts{v}) ? "INFO" : "WARNING";
my $Logger = HUBDesign::Logger->new(level => $opts{v});

#Handle Simple Numeric Options 
$opts{e} = ProcessNumericOption($opts{e},$DEFAULT{e},undef,undef,0,"e-value");
$opts{l} = ProcessNumericOption($opts{l},$DEFAULT{l},1,undef,1,"Probe Length");
$opts{t} = ProcessNumericOption($opts{t},$DEFAULT{t},0,undef,0,"Max Threads");

#Ensure Threads is within system limits
$opts{t} = ValidateThreadCount($opts{t});

#Handle Command line arguments
my %FileDict= (BRInfo => shift @ARGV);
my @BlastDBList = @ARGV;

$Logger->Log("Bait Region File ($FileDict{BRInfo}) doe not exits","ERROR") unless(-e$FileDict{BRInfo});

#Prep parameters for logging
my %params;
@params{('Temp Directory', 'Probe Length', 'Threads', 'Input Data')} = (@opts{qw(d l t)}, $FileDict{BRInfo});
$params{'Keep Temp Files'} = ($opts{k}) ? 'TRUE' : 'FALSE';
$params{'BlastDB Count'} = scalar(@BlastDBList);

#============================================================
#Main

    $Logger->LogParameters(%params);
    my %BaitRegionDict = LoadBaitRegions($FileDict{BRInfo});
    $FileDict{Fasta_CBR} = OutputFasta(\%BaitRegionDict,"CBR");
    unless($opts{f} eq 'no'){
        my %ExIVDict = RunLCRFilter($FileDict{Fasta_CBR});
        ExcludeIntervals(\%BaitRegionDict,\%ExIVDict);
    }
    for(my $i = 0; $i < @BlastDBList; $i++){
        last unless(scalar (keys %BaitRegionDict));
        my $blastDB = $BlastDBList[$i];
        my $label= "FBR$i";
        $FileDict{"Fasta_$label"} = OutputFasta(\%BaitRegionDict,$label);
        my %ExIVDict = RunBlastFilter($FileDict{"Fasta_$label"},$blastDB);
        ExcludeIntervals(\%BaitRegionDict,\%ExIVDict);
    }
   OutputResults([values %BaitRegionDict]); 

#============================================================
#Subroutines

#Parses A Candidate Bait Region File loading those regions into memory
#Inputs - Path to a tab delimited file with 5 columns: bait region id, taxon id, cluster id, position
   #within cluster, and the regions sequence
#Output - A hash with keys of bait region id and values of HUBDesign::BaitRegion objects
sub LoadBaitRegions($){
    $Logger->Log("Loading candidate bait regions...","INFO");
    my $fh = OpenFileHandle(shift,"Candidates","ERROR");
    my %BRDict;
    my %TaxonSet;
    while(my $line = <$fh>){
        chomp($line);
        my ($br_id, $taxon_id, $clust_id, $pos, $seq) = split(/\t/,$line);
        if(exists $BRDict{$br_id}){
            $Logger->Log("Ignoring duplicate entry for bait region $br_id on line $.","WARNING");
            next;
        }
        $TaxonSet{$taxon_id} = 1;
        $BRDict{$br_id} = BaitRegion->new(taxon_id => $taxon_id, clust_id => $clust_id, pos => $pos, seq => $seq);
    }
    my $message = sprintf("Loaded %d candidate bait regions for %d taxa",scalar(keys %BRDict),scalar(keys %TaxonSet));
    $Logger->Log($message,"INFO");
    return %BRDict;
}

#Given a dictionary of candidate bait regions, outputs a fasta formated file
#Inputs - reference to a hash of HUBDesign::BaitRegion objects
#       - a label for the fasta file, mostly for logging purposes
#Output - A File::Temp object where the fasta sequence is written
sub OutputFasta($$){
    my ($BRDictRef, $label) = @_;
    $Logger->Log("Writing candidate bait regions as fasta...","INFO");
    my $tmpFile = File::Temp->new(DIR => $opts{d}, TEMPLATE => "HUBDesign_${label}_XXXX", suffix => ".fna",UNLINK => ($opts{k} + 1) % 2);
    my $SeqOObj = Bio::SeqIO->new(-file => ">$tmpFile", -format => 'fasta');
    my $written = 0;
    my $length = 0;
    foreach my $br_id (keys %{$BRDictRef}){
        my $seqObj = Bio::Seq->new(-id => $br_id, -seq => $BRDictRef->{$br_id}->seq);
        $SeqOObj->write_seq($seqObj);
        $written++;
        $length += $seqObj->length;
    }
    $Logger->Log("Wrote $written fasta entries totalling $length bp to $tmpFile","INFO");
    return $tmpFile;

}

#Given A fasta formatted set of sequences runs the lcr filter to identify lcr containing intervals
#Inputs - A path/File::Temp object to a fasta formattted file
#Output - A hash with keys of bait region id and values two-element array refs [start, end]
sub RunLCRFilter($){
    my $file = shift;
    $Logger->Log("Filtering candidate bait regions for low complexity sequences...","INFO");
    my $cmd = "$_LCR_EXECUTABLE $opts{f} $file";
    my $fh = OpenFileHandle("$cmd |","LCR cmd","WARNING");
    unless($fh){
        $Logger->Log("LCR Filtration Skipped","INFO");
        return ();
    }
    my %ExIVDict;
    my $count = 0;;
    my $length = 0;
    while(my $line = <$fh>){
        chomp($line);
        my ($br_id, $start, $end) = split(/\t/,$line);
        #sdust is 0-indexed, desire 1-indexed
        $start++;
        $end++;
        $ExIVDict{$br_id} = [] unless(exists $ExIVDict{$br_id});
        push(@{$ExIVDict{$br_id}},[$start, $end]);
        $count++;
        $length += $end - $start + 1;
    }
    close($fh);
    $Logger->Log("Identified $count lcrs totalling $length bp","INFO");
    return %ExIVDict;
}

#Given a hash of excludion intervals, modifies a hash of bait regions, breaking up larger regions which
#contain intervals
#Inputs - A reference to a hash of HUBDesign::BaitRegion objects
#       - A reference to a hash of [start, end] array refs
#Output - None, the bait region hash is modified in place
sub ExcludeIntervals($$){
    my ($BRDictRef,$ExIVListRef) = @_;
    foreach my $br_id (keys %{$ExIVListRef}){
        my @new_baitReg = $BRDictRef->{$br_id}->exclude($ExIVListRef->{$br_id});
        delete($BRDictRef->{$br_id});
        for(my $i = 0; $i < @new_baitReg; $i++){
            if(length($new_baitReg[$i]->seq) >= $opts{l}){
                $BRDictRef->{$br_id."_".$i} = $new_baitReg[$i];
            }
        }
    }
    my %TaxonSet;
    foreach my $brObj (values %{$BRDictRef}){
        $TaxonSet{$brObj->taxon_id} = 1;
    }
    my $message = sprintf("There are %d bait regions from %d taxa remaining",scalar(keys %{$BRDictRef}),scalar(keys %TaxonSet));
    $Logger->Log($message,"INFO");
}

#Given a fasta formatted file and a database, runs blast to identify intervals for exclusion
#Inputs - A path/File::Temp object to a fasta formattted file
#       - A path to a blast database
#Output - A hash with keys of bait region id and values two-element array refs [start, end]
sub RunBlastFilter($$){
    my ($fasta,$db) = @_;
    $Logger->Log("Filtering candidate bait regions against $db ...","INFO");
    my $addparam = ($opts{x} eq 'none') ? '' : $opts{x};
    my $cmd = "$_BLAST_EXECUTABLE -task blastn $addparam -db $db -num_threads $opts{t} -query $fasta";
    $cmd .= " -evalue $opts{e} -outfmt '6 qseqid length pident qstart qend' -dust 'no'";
    $Logger->Log("BlastCmd: $cmd","INFO");
    my $fh = OpenFileHandle("$cmd |","Blast cmd","ERROR");
    my %ExIVDict;
    while(my $line = <$fh>){
        chomp($line);
        my ($br_id,$length,$pident, $start, $end) = split(/\t/,$line);
        next unless ($length >= $opts{m});
        next unless ($pident >= $opts{i});
        $ExIVDict{$br_id} = [] unless(exists $ExIVDict{$br_id});
        ($start,$end) = ($end,$start) if($start > $end);
        push(@{$ExIVDict{$br_id}},[$start, $end]);
    }
    close($fh);
    while (my ($br_id,$ivListRef) = each %ExIVDict){
        @{$ivListRef} = sort {$a->[0] <=> $b->[0]} @{$ivListRef};
        
        for(my $i = 1; $i < @{$ivListRef}; $i++){
            if($ivListRef->[$i-1]->[1] >= $ivListRef->[$i]->[0]){
                my $newend = $ivListRef->[$i-1]->[1];
                $newend =  $ivListRef->[$i]->[1] if($ivListRef->[$i]->[1] > $newend);
                $ivListRef->[$i-1]->[1] = $newend;
                splice(@{$ivListRef},$i,1);
                $i--;
            } 
        }
    }
    my $count = 0;;
    my $length = 0;
    while (my ($br_id,$ivListRef) = each %ExIVDict){
        foreach my $ivRef (@{$ivListRef}){
            $count++;
            $length += $ivRef->[1] - $ivRef->[0] + 1;
        }
    }
    $Logger->Log("Identified $count intervals totalling $length bp\n","INFO");
    return %ExIVDict;

}

#Given a set of bait regions, outputs them
#Inputs - A reference to a hash with keys of bait region id and values of HUBDesign::BaitRegion objects
#Output - None, prints to STDOUT
sub OutputResults($){
    my $BRListRef = shift;
    $Logger->Log("Outputting Bait Regions...","INFO");
    for(my $i = 0; $i < @{$BRListRef}; $i++){
        my $reg_id = sprintf("FBR_%0*d",$_ID_PADDING,$i);
        print  "$reg_id\t",$BRListRef->[$i]->toStr,"\n";
    }
    $Logger->Log("Done","INFO");
}
