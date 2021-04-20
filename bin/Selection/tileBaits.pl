#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Scalar::Util qw(looks_like_number);
use Getopt::Std;
use Class::Struct;
use POSIX;
use sort 'stable';

sub ProcessNumericOption($$$$$$);
sub ProcessDelimStr_ListFile($$$$);
sub PrintREGION($$);
sub GetSortedRegions();
sub sortMaxSpread(@);
sub MinDistance($@);
sub GetTilingDens($);                       #Usage: my TilingObj = GetTiling(seqLen);
sub CalculateSpacing($$);
sub CalculateNumBaits($$;$);
sub GetBaitCount($@);                       #Usage sub GetMaximumBaitCount(('Min','Max','Tgt',$Int),@ListOfOrderedRegIDListRefs);
sub GetRegionDensityByMode($$);             #Usage my $density = GetRegionDensityByMode($Mode,$regID);
sub ApplyMaxRegionFilter(\%;$);
sub OutputBaits(@);                         #Usage OutputBaits(@ListOfOrderedRegIDListRefs)
sub DetermineTilingClass();                 #Usage my @OrgTilingClass = DetermineTilingClass();
sub GetBaitCountsPerOrg();
sub SetRegionTargetTiling(%);               #Usage SetRegionTargetTiling(%OrgTilingClass);
sub CalculateWeightofRegionWithinOrg($$$);  #Usage: fncname($IdealBaitsPerOrg,$orgID,$regID);
sub CalculateTilingDensity($$);             #Usage: fncname($regLen,$numBaits);
sub TrimBaits($$);                           #Usage: $BaitsTrimmed = TrimBaits($NumBaitstoTrim);
sub TrimBaitsFromOrganisms($$%);             #Usage: $BaitsTrimmed = TrimBaitsFromOrganism($densityMode,(Org1 => NumToRemove))
sub ApplyMaxBaitPerOrgFilter(;$);           #Usage: 
sub OptimizeMaxBaitsPerOrg();               
sub OptimizeMaxRegionsPerOrg();

sub PrintUpdate($$);


#min = Minimum Possible Tiling Density
#max = Maximum Possible Tiling Density
#tgt = Target Density
struct(TILING => {min => '$', max => '$', tgt => '$'});

#Org = An array of all organisms this bait region applies to
#Tgt = The type of sequence this region targets
#Seq = The actual sequence of the bait region
#Pos = The position within the target sequence this region begins at
#TDn = A TILING structure
#Cls = The tiling class of this region, one of Min, Max, or Med, the last can be progressively tiled
#Idx = The position within the region to begin tiling at
#Len = The length of the tiling region of the bait region
#Ex. A 400bp Sequence which has been trimmed to 4 baits at 4x coverage would have: Len=170, Idx=115
#Ex. A 150bp Sequence which has been trimmed to 37 baits at 35x coverage would have Len=149, Idx=0
struct(REGION => {Org => '@', Tgt => '$', Seq => '$', Pos => '$',
        TDn => 'TILING', Cls => '$', Idx => '$', Len => '$'});

my $Command = join(' ',($0,@ARGV));

my %opts;
getopts('DOb:d:l:n:p:r:aPSvh',\%opts);

if(exists $opts{h} or @ARGV < 2){
    my $Usage = "@{[basename($0)]} [-options] RegionInfo RegionSeq ... > TiledBaits.fna";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of Bait Regions and some metadata, outputs a set of baits which\n".
            "\t\toptimally tile the regions.\n".
            "\tOptimal is baits centered within regions, with maximally balanced baits per species\n".
            "\t\twith maximally spaced small regions given priority\n".
            "===Usage\n".
            "\t$Usage\n".
            "===Required Inputs\n".
            "\tRegionInfo\tA tab delimited file with 4 columns:\n".
            "\t\tRegionID (matches fasta headers), Target Organism, Target Seq within organism , Position within Seq\n".
           "\tRegionSeq\tA fasta formated file or set of files\n".
           "===Options\n".
           "\t-b INT[0,∞)\tThe maximum number of baits for any target organism [Default: 0 -> None]\n".
           "\t-d INT[1,l]\tThe minimum tiling density for any bait region [Default: 1x]\n". 
           "\t-l INT[1,∞)\tThe length of baits to output [Default: 75bp]\n".
           "\t-n INT[0,∞)\tThe maximum number of baits in the output [Default: 0 -> None]\n".
           "\t-p string\tA comma separated string or list file detailing the order of prescedence\n".
           "\t\tfor target sequences. By default all are considered equal\n".
           "\t-r INT[0,∞)\tThe maximum number of regions for any target organism [Default: 0 -> None]\n".
           "===Flags\n".
           "\t-a\tAuto-detect Maximum and Regions Baits per Organism to improve Proportional Tiling\n".
           "\t-P\tDisable Proportional tiling, and instead find a maximum constant tiling\n".
           "\t-S\tSkip Region Sorting\n".
           "\t-v\tVerbose Output\n".
           "\t-h\tDisplay this message and exit\n";
    } else {
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

my ($MetadataFile,@FastaFileList) =  @ARGV;

##INTERNAL VARS

my $bDebug = (exists $opts{D}) ? 1 : 0;
my $OutMode = ('TaxaInfo','Baits')[exists $opts{O} ? 0 : 1];
my $minBaitsPerOrg = 30;

## Process Input
my $maxBaitsPerOrg = ProcessNumericOption($opts{b},0,0,undef,1,"Maximum Baits per Organism");
my $baitLen = ProcessNumericOption($opts{l},75,1,undef,1,"BaitLength");
my $minDensity = ProcessNumericOption($opts{d},1,1,$baitLen,1,"Minimum Tiling Density");
my $maxBaits = ProcessNumericOption($opts{n},0,0,undef,1,"Maximum Baits");
my $maxRegPerOrg = ProcessNumericOption($opts{r},0,0,undef,1,"Maximum Regions per Organism");
my @TargetSeqOrder = ProcessDelimStr_ListFile($opts{p},",",[],"Target Sequence Order");
if(exists $opts{a}){
    if(exists $opts{P}){
        warn "[WARNING] Auto-detect maxBaitsPerOrg option is useless if proportional tiling is disabled\n";
    }
    if($maxBaitsPerOrg){
        warn "[WARNING] Auto-detect maxBaitsPerOrg overules user specification: Auto-detecting\n";
        $maxBaitsPerOrg = 0;
    }
    if($maxRegPerOrg){
        warn "[WARNING] Auto-detect maxRegPerOrg overules user specification: Auto-detecting\n";
        $maxRegPerOrg = 0;
    }
}
my $bRegSort = exists $opts{S} ? 0 : 1;

my $InitTime = time;
my $LastTime = $InitTime;
my ($sec,$min,$hour,$mday,$mon,$year) = localtime($InitTime);
my @abbMon = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
$year += 1900;
if(exists $opts{v}){
    print STDERR "[INFO] Initialized at $hour:$min:$sec on $mday $abbMon[$mon] $year\n"; 
    my @tmp = ("=") x 80;
    print STDERR join("",@tmp),"\n";
    my $Str = "[INFO] Running Tiling Regions with the following parameters:\n";
    $Str .= "\tcmd:\t\t$Command\n";
    $Str .= "\tmaxBaitsPerOrg:\t$maxBaitsPerOrg\n" if($maxBaitsPerOrg);
    $Str .= "\tmaxBaitsPerOrg:\tAuto-detect\n" if(exists $opts{a});
    $Str .= "\tbaitLen:\t$baitLen\n";
    $Str .= "\tminDensity:\t$minDensity\n";
    $Str .= "\tmaxBaits:\t$maxBaits\n" if($maxBaits);
    $Str .= "\tmaxRegPerOrg:\t$maxRegPerOrg\n" if($maxRegPerOrg);
    $Str .= "\tmaxRegPerOrg:\tAuto-detect\n" if(exists $opts{a});
    $Str .= "\tTargetSeqOrder:\t@{[join(',',@TargetSeqOrder)]}\n" if(@TargetSeqOrder);
    $Str .= "\tTilingMode:\t@{[(exists $opts{P} ? 'Maximal' : 'Proportional')]}\n";
    $Str .= "\tMetadata:\t$MetadataFile\n";
    $Str .= "\tSequenceFiles:\t@{[join(',',@FastaFileList)]}\n";
    print STDERR $Str;
    PrintUpdate("Loading Region Metadata",0);
}


## Load Metadata
my %BaitRegionInfo; # Hash with keys of region ids and values of REGION
if(open(my $fh,$MetadataFile)){
    while(my $line = <$fh>){
        chomp($line);
        my ($id,$orgstr,$target,$pos) = split(/\t/,$line);
        if(exists $BaitRegionInfo{$id}){
            warn "[WARNING] Multiple Metadata entries present for $id: Only the first is used\n";
        } else{
            $BaitRegionInfo{$id} = REGION->new(Org => [split(",",$orgstr)],Tgt => $target,Pos=>$pos,Idx=>0);
        }
    }
    close($fh);
} else {
    die "[ERROR] Could not open Metadata file - $MetadataFile:$!\n";
}

if(exists $opts{v}){
    PrintUpdate("Loaded Metadata for @{[scalar(keys %BaitRegionInfo)]} Bait Regions",1);
    PrintUpdate("Loading Sequences",0);
}

## Load Sequences
my $seqCount =0;
foreach my $file (@FastaFileList){
    die "ERROR: Fasta File - $file does not exist\n" unless( -e $file);
    my $IN = Bio::SeqIO->new(-file => $file, -format => 'fasta');
    while (my $seqObj = $IN->next_seq){
        my $id = $seqObj->id;
        if(exists $BaitRegionInfo{$id}){
            $seqCount++;
            $BaitRegionInfo{$id}->Seq($seqObj->seq);
            $BaitRegionInfo{$id}->Len(length($seqObj->seq));
            $BaitRegionInfo{$id}->TDn(GetTilingDens(length($seqObj->seq)));
        } else {
            warn "[WARNING] No metadata provided for sequence $id: Ignoring this sequence\n";
        }
    }
}

##Ignore any metadata from regions with no sequence
foreach my $regID (keys %BaitRegionInfo){
    delete $BaitRegionInfo{$regID} unless (defined $BaitRegionInfo{$regID}->Len);
}

if(exists $opts{v}){
    PrintUpdate("Loaded $seqCount Sequences",1);
    PrintUpdate("Sorting Regions within Organisms",0);
}

## MAIN PROCESSING
my $CalcTime = time;

#Define order of regions within taxa
my $SortTime = time;
my %RegionOrder = GetSortedRegions();

$LastTime = $SortTime;
PrintUpdate("Regions Sorted",1) if(exists $opts{v});
if($maxRegPerOrg){
    PrintUpdate("Applying Maximum Region Per Organism Filter",0) if(exists $opts{v});
    ApplyMaxRegionFilter(%RegionOrder);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach((values %RegionOrder));
    if(exists $opts{v}){
        PrintUpdate("Filtered out @{[scalar(keys %BaitRegionInfo) - scalar(keys %UniqRegID)]} Regions",1);
    }
}
my $AbsMaxBaits = GetBaitCount('Max',(values %RegionOrder));
if($maxBaitsPerOrg){
    PrintUpdate("Applying Maximum Bait Per Organism Filter",0) if(exists $opts{v});
    my $PreMaxBaits = $AbsMaxBaits;
    ApplyMaxBaitPerOrgFilter();
    $AbsMaxBaits = GetBaitCount('Max',(values %RegionOrder));
    PrintUpdate("Filtered out @{[$PreMaxBaits - $AbsMaxBaits]} Baits",1) if(exists $opts{v});
}
my %OrgTilingClass;
if($maxBaits && $AbsMaxBaits > $maxBaits){
    if(exists $opts{P}){ #Select maximal tiling keeping # Baits below thresh
        PrintUpdate("Determining Maximal Constant Tiling Density",0) if(exists $opts{v});
        my $density = $minDensity;
        $density++ while(GetBaitCount($density+1,(values %RegionOrder)) <= $maxBaits);
        foreach my $regID (keys %BaitRegionInfo){
            my $min = $BaitRegionInfo{$regID}->TDn->min;
            my $max = $BaitRegionInfo{$regID}->TDn->max;
            my $tgt = ($density < $min) ? $min : $density;
            $tgt = ($tgt > $max) ? $max : $tgt;
            $BaitRegionInfo{$regID}->TDn->tgt($tgt);
        }
        PrintUpdate("Found Maximal Constant Tiling Density at ${density}x",1) if(exists $opts{v});
        (undef,%OrgTilingClass) = DetermineTilingClass();
    } else { #Scale Tiling to balance #Baits per Org
         if(exists $opts{a}){
             PrintUpdate("Optimizing Max Regions Per Organism",0) if(exists $opts{v});
             my $OptimTime = time;
             OptimizeMaxRegionsPerOrg();
             $LastTime = $OptimTime;
             PrintUpdate("Found Optimum at $maxRegPerOrg Regions Per Organism",1) if(exists $opts{v});
             PrintUpdate("Optimizing Max Baits Per Organism",0) if(exists $opts{v});
             $OptimTime = time;
             OptimizeMaxBaitsPerOrg();
             $LastTime = $OptimTime;
             PrintUpdate("Found Optimum at $maxBaitsPerOrg Baits Per Organism",1) if(exists $opts{v});
         }
         PrintUpdate("Selecting Proportional Tiling Densities",0) if(exists $opts{v});
         my $Ideal;
         ($Ideal,%OrgTilingClass) = DetermineTilingClass();
         SetRegionTargetTiling(%OrgTilingClass);
         if(exists $opts{v}){
             $Ideal = sprintf("%.0f",$Ideal);
             PrintUpdate("Balanced Tiling Densities Target $Ideal Baits per Organism",0);
         }
    }
    #foreach (values %BaitRegionInfo){
    #    PrintREGION($_,*STDERR);
    #}
    my $curBaits = GetBaitCount('Tgt',(values %RegionOrder));
    while($curBaits > $maxBaits){
        PrintUpdate("Currently have $curBaits: Trimming Baits down to $maxBaits",0) if(exists $opts{v});
        my %Before = GetBaitCountsPerOrg();
        TrimBaits('Tgt',$curBaits - $maxBaits);
        my %After = GetBaitCountsPerOrg();
        for my $org (keys %Before){
            foreach my $cls (('min','max','tgt')){
                if($Before{$org}->{$cls} < $After{$org}->{$cls}){
                    print STDERR "$org:$cls: $Before{$org}->{$cls} -> $After{$org}->{$cls}\n"
                }
            }
        }
        my $nowBaits = GetBaitCount('Tgt',(values %RegionOrder));
        if(exists $opts{v}){
            PrintUpdate("Trimmed @{[$curBaits - $nowBaits]} Baits",1);
        }
        $curBaits = $nowBaits;
    }
} elsif($maxBaits) {
    warn "[INFO] Maximum Possible Baits is ≤ specified maximum : All baits to be output at maximum density\n";
    foreach my $RegObj (values %BaitRegionInfo){
        $RegObj->TDn->tgt($RegObj->TDn->max);
    }
}
$LastTime=$CalcTime;
if(exists $opts{v}){
    my $BaitCount = GetBaitCount('Tgt',(values %RegionOrder));
    PrintUpdate("Calculations Complete",1);
    PrintUpdate("Outputting $BaitCount Baits",0);
}
if($OutMode eq 'Baits'){
    OutputBaits((values %RegionOrder));
} elsif($OutMode eq 'TaxaInfo') {
    my %BaitCountPerOrg = GetBaitCountsPerOrg();
    foreach my $org (sort {$BaitCountPerOrg{$a}->{tgt} <=> $BaitCountPerOrg{$b}->{tgt}} keys %BaitCountPerOrg){
        print "$org\t$BaitCountPerOrg{$org}->{tgt}\t$OrgTilingClass{$org}\n";
        #if($org eq '478820'){
        #    foreach my $regID (@{$RegionOrder{$org}}){
        #        print STDERR"ID: $regID\n";
        #        PrintREGION($BaitRegionInfo{$regID},*STDERR);
        #    }
        #}
    }
}
if(exists $opts{v}){
    PrintUpdate("Done Outputing Baits",1);
    $LastTime=$InitTime;
    PrintUpdate("Processing complete",1);
    my @tmp = ("=") x 80;
    print STDERR join("",@tmp),"\n";
}


## FUNCTION DEFINITIONS

sub ProcessNumericOption($$$$$$){
    print STDERR "[DEBUG] ProcessNumericOption\n" if($bDebug);
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


sub ProcessDelimStr_ListFile($$$$){
    print STDERR "[DEBUG] ProcessDelimStr_ListFile\n" if($bDebug);
    my($val,$delim,$default,$varName) = @_;
    return @{$default} unless(defined $val);
    my @valList = split($delim,$val);
    my @Output;
    foreach my $file (@valList){
        if(-e $file){
            open(my $fh, $file) or die "[ERROR] Could not open $varName file - $file: $!\n"; 
            my @lines = <$fh>;
            chomp(@lines);
            push(@Output,@lines);
            close($fh);
        } else {
            push(@Output,$file);
        }
    }
    return @Output;
}


sub PrintREGION($$){
    print STDERR "[DEBUG] PrintREGION\n" if($bDebug);
    my ($Reg,$fh) = (@_);
    my $oldfh = select($fh);
    print "Org: ".join(',',@{$Reg->Org})."\n";
    print "Tgt: ".$Reg->Tgt."\n";
    print "Pos: ".$Reg->Pos."\n";
    if(length($Reg->Seq) > 20){
        print "Seq: ".substr($Reg->Seq,0,10)." ... ".(length($Reg->Seq) - 20).
            " Residues  ... ".substr($Reg->Seq,-10)."\n";
    } else {
        print "Seq: ".$Reg->Seq."\n";
    }
    print "Len: ".$Reg->Len."\n";
    print "Idx: ".$Reg->Idx."\n";
    print "Cls: ".$Reg->Cls."\n" if(defined $Reg->Cls);
    print "TDn:\n";
    print "\tmax: ".$Reg->TDn->max."\n" if(defined $Reg->TDn->max);
    print "\tmin: ".$Reg->TDn->min."\n" if(defined $Reg->TDn->min);
    print "\ttgt: ".$Reg->TDn->tgt."\n" if(defined $Reg->TDn->tgt);
    select($oldfh);
}

sub GetSortedRegions() {
    print STDERR "[DEBUG] GetSortedRegions\n" if($bDebug);
    my %RegionOrder;
    ##Assign each region to its organism
    foreach my $regID (sort keys %BaitRegionInfo){
        my $RegionObj = $BaitRegionInfo{$regID};
        foreach my $org (@{$RegionObj->Org}){
            unless(exists $RegionOrder{$org}){
                $RegionOrder{$org} = [];
            }
            push(@{$RegionOrder{$org}},$regID);
        }
    }
    return %RegionOrder unless($bRegSort);
    my %TargetSeqRank;
    if(@TargetSeqOrder){
        for(my $i = 0; $i < @TargetSeqOrder; $i++){
            $TargetSeqRank{$TargetSeqOrder[$i]} = $i;
        }
    }
    my $OrgCounter = 0;
    my $OrgCount = scalar(keys %RegionOrder);
    while (my ($org,$idListRef) = each %RegionOrder){
        my $elapsed = time - $LastTime;
        $OrgCounter++;
        if($opts{v} and $elapsed > 10){
            PrintUpdate("Now sorting $org ($OrgCounter/$OrgCount)",1);
        }
        ##Sort first by Size - Prioritizing Small Regions
        #PrintUpdate("Size Sort $org ($OrgCounter/$OrgCount)",1) if(exists $opts{v});
        @{$idListRef} = sort {$BaitRegionInfo{$a}->Len <=> $BaitRegionInfo{$b}->Len} @{$idListRef};
        ##Sort such that maximal spacing is maintained (Stable sort: Size to break ties)
        #PrintUpdate("Position Sort $org ($OrgCounter/$OrgCount)",1) if(exists $opts{v});
        my @idxPosList;
        for(my $i = 0; $i < @{$idListRef}; $i++){
            push(@idxPosList,{idx => $i, pos => $BaitRegionInfo{$idListRef->[$i]}});
        }
        @idxPosList = sortMaxSpread(@idxPosList);
        @{$idListRef} = (@{$idListRef})[map {$_->{idx}} @idxPosList];
        ###Sort in order of target precedence (Stable sort : Spacing then Size to break ties)
        #PrintUpdate("Target Sort $org ($OrgCounter/$OrgCount)",1) if(exists $opts{v});
        if(@TargetSeqOrder){
            foreach  my $regID (@{$idListRef}){
                unless(exists $TargetSeqRank{$BaitRegionInfo{$regID}->Tgt}){
                    $TargetSeqRank{$BaitRegionInfo{$regID}->Tgt} = scalar(@TargetSeqOrder) + 1;
                }
            }
            @{$idListRef} = sort {$TargetSeqRank{$BaitRegionInfo{$a}->Tgt} <=> $TargetSeqRank{$BaitRegionInfo{$b}->Tgt}} @{$idListRef}
        }; 
        ##Sort in order of specificiy - Prioritizing fewer target organisms (Stable Sort)
        #PrintUpdate("Target Sort $org ($OrgCounter/$OrgCount)",1) if(exists $opts{v});
        @{$idListRef} = sort {scalar(@{$BaitRegionInfo{$a}->Org}) <=> scalar(@{$BaitRegionInfo{$b}->Org})} @{$idListRef};
    }
    return %RegionOrder;
}

sub sortMaxSpread(@){
    print STDERR "[DEBUG] sortMaxSpread\n" if($bDebug);
    my @regList = @_;
    @regList = sort {$a->{pos} <=> $b->{pos}} @regList;
    my @pairList = ([0,$#regList]);
    my @output = (@regList[@{$pairList[0]}]);
    while (@pairList) {
        my $pair = shift(@pairList);
        my $iRange = $pair->[1]-$pair->[0];
        if($iRange <= 2){
            if($iRange == 2){
                push(@output,$regList[$pair->[0]+1]);
            }
            next;
        }
        my $range = $regList[$pair->[1]]->{pos} - $regList[$pair->[0]]->{pos};
        my $m = $pair->[0]+1;
        my $dist = $regList[$m]->{pos} - $regList[$pair->[0]]->{pos};
        until($dist >= $range / 2){
            $m++;
            $dist += $regList[$m]->{pos} - $regList[$m -1]->{pos};
        }
        $m-- if($m == $pair->[1]);
        push(@output,$regList[$m]);
        push(@pairList,[$pair->[0],$m]);
        push(@pairList,[$m,$pair->[1]]);
    }
    return @output;
}

#sub sortMaxSpread(@){
#    print STDERR "[DEBUG] sortMaxSpread\n" if($bDebug);
#    my @regList = @_;
#    @regList = sort {$a->{pos} <=> $b->{pos}} @regList;
#    my $m = int($#regList/2);
#    my @posList = ($regList[$m]->{pos});
#    my @output = ($regList[$m]);
#    if($m > 0 and $m < $#regList){
#        @regList = @regList[(0 .. ($m -1), ($m + 1) .. $#regList)];
#    } elsif($m == 0 and @regList == 1){
#        shift(@regList);
#    } elsif($m == 0){
#        @regList = @regList[1 .. $#regList];
#    } else{
#        @regList = @regList[0 .. ($#regList -1)];
#    }
#    while(@regList){
#        @regList = sort {MinDistance($b->{pos},@posList) <=> MinDistance($a->{pos},@posList)} @regList;
#        push(@posList,$regList[0]->{pos});
#        push(@output,shift(@regList));
#    }
#    return @output;
#}


#sub MinDistance($@){
#    print STDERR "[DEBUG] MinDistance\n" if($bDebug);
#    my($pos,@posList) = @_;
#    my $min = undef;
#    foreach my $loc (@posList){
#        $min = abs($pos - $loc) if(!defined($min) or abs($pos - $loc) < $min);
#    }
#    return $min;
#}

sub GetTilingDens($){
    print STDERR "[DEBUG] GetTilingDens\n" if($bDebug);
    my ($regLen) = @_;
    my $maxDen = ($regLen >= 2*$baitLen-1) ? $baitLen : $regLen-$baitLen+1;
    my $minDen = ($maxDen < $minDensity) ? $maxDen : $minDensity;
    return TILING->new(max => $maxDen, min => $minDen);
}

sub CalculateSpacing($$){
    print STDERR "[DEBUG] CalculateSpacing(@_) @ @{[ caller ]}\n" if($bDebug);
    my ($regLen,$density) = @_;
    return $baitLen if ($density == 1);
    return ($regLen < 2*$baitLen) ? floor(($regLen - $baitLen)/($density-1)) : floor(($baitLen - 1)/($density - 1));

}

sub CalculateNumBaits($$;$){
    print STDERR "[DEBUG] CalculateNumBaits(@_) @ @{[ caller ]}\n" if($bDebug);
    my ($regLen,$density,$regID) = @_;
    my $spacing = CalculateSpacing($regLen,$density);
    if(!$spacing && $bDebug){
        print STDERR "ID: $regID\n";
        PrintREGION($BaitRegionInfo{$regID},*STDERR);
    }
    return ceil(($regLen - $baitLen + 1) / $spacing);
}

sub GetBaitCount($@){
    print STDERR "[DEBUG] GetBaitCount(@_) @ @{[ caller ]}\n" if($bDebug);
    my $densMode = shift(@_);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(@_);
    my $count = 0;
    foreach my $regID (keys %UniqRegID){
        my $density = GetRegionDensityByMode($densMode,$regID); 
        $count += CalculateNumBaits($BaitRegionInfo{$regID}->Len,$density,$regID);
    }
    return $count;
}

sub GetBaitCountsPerOrg(){
    print STDERR "[DEBUG] GetBaitCountsPerOrg\n" if($bDebug);
    my %BaitCountsPerOrg;
    while( my ($org,$regionOrderRef) = each %RegionOrder){
        my %Count = (min => 0, max => 0, tgt => 0);
        foreach my $regID (@{$regionOrderRef}){
            $Count{min} += CalculateNumBaits($BaitRegionInfo{$regID}->Len,$BaitRegionInfo{$regID}->TDn->min);
            $Count{max} += CalculateNumBaits($BaitRegionInfo{$regID}->Len,$BaitRegionInfo{$regID}->TDn->max);
            if(defined $BaitRegionInfo{$regID}->TDn->tgt){
                $Count{tgt} += CalculateNumBaits($BaitRegionInfo{$regID}->Len,$BaitRegionInfo{$regID}->TDn->tgt);
            }
        }
        $BaitCountsPerOrg{$org} = \%Count;
    }
    return %BaitCountsPerOrg;
}

sub GetRegionDensityByMode($$){
    print STDERR "[DEBUG] GetRegionDensityByMode(@_) @{[ caller ]}\n" if($bDebug);
    my($densMode,$regID) = @_;
    my $density;
    if($densMode eq 'Max'){
        $density = $BaitRegionInfo{$regID}->TDn->max;
    } elsif($densMode eq 'Min'){
        $density = $BaitRegionInfo{$regID}->TDn->min;
    } elsif($densMode eq 'Tgt'){
        $density = $BaitRegionInfo{$regID}->TDn->tgt;
    } elsif(looks_like_number($densMode)){
        $density = $densMode;
        if($density < $BaitRegionInfo{$regID}->TDn->min){
            $density = $BaitRegionInfo{$regID}->TDn->min;
        }elsif($density > $BaitRegionInfo{$regID}->TDn->max){
            $density = $BaitRegionInfo{$regID}->TDn->max;
        }
    } else {
        die "Unknown density Mode: $densMode\n";
    }
    return $density;
}

sub ApplyMaxRegionFilter(\%;$){
    print STDERR "[DEBUG] ApplyMaxRegionFilter\n" if($bDebug);
    #maxRegPerOrg is assumed to be non-zero
    my($RegionOrderRef,$bQuiet) = @_;
    $bQuiet = 0 unless(defined $bQuiet);
    #ID All Organisms with > Threshold Regions
    my %OrgAboveThresh;
    while(my ($org,$orderRef) = each %{$RegionOrderRef}){
        if(@{$orderRef} > $maxRegPerOrg){
            $OrgAboveThresh{$org} = @{$orderRef} - $maxRegPerOrg;
        }
    }
    foreach my $org (sort {$OrgAboveThresh{$b} <=> $OrgAboveThresh{$a}} keys %OrgAboveThresh){
        for(my $i = $#{$RegionOrderRef->{$org}}; $i >= 0 and exists $OrgAboveThresh{$org}; $i--){
            my $regID = $RegionOrderRef->{$org}->[$i];
            my @OrgList = @{$BaitRegionInfo{$regID}->Org};
            next if(scalar(@OrgList) == 1 and $i == 0); #Do not remove the last organism specific region
            my $bRemove = 1;
            for(my $j = 0; $j < @OrgList and $bRemove; $j++){
                $bRemove = 0 unless(exists $OrgAboveThresh{$OrgList[$j]});
            }
            if($bRemove){
                foreach my $tgtorg (@OrgList){
                    for (my $j = 0; $j < @{$RegionOrderRef->{$tgtorg}}; $j++){
                        if($RegionOrderRef->{$tgtorg}->[$j] eq $regID){
                            splice(@{$RegionOrderRef->{$tgtorg}},$j,1);
                            last;
                        }
                    }
                    unless(--$OrgAboveThresh{$tgtorg}){
                        delete $OrgAboveThresh{$tgtorg};
                    }
                }
            }
        }
        if(exists $OrgAboveThresh{$org} and !$bQuiet){
            warn "[WARNING] Could not acheive at most $maxRegPerOrg regions for $org, due to region sharing: acheived @{[$maxRegPerOrg + $OrgAboveThresh{$org}]} regions\n";
        }
    }
}


sub OutputBaits(@){
    print STDERR "[DEBUG] OutputBaits\n" if($bDebug);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(@_);
    my $OUT = Bio::SeqIO->newFh(-format => "fasta");
    foreach my $regID (keys %UniqRegID){
        my $density = (defined $BaitRegionInfo{$regID}->TDn->tgt) ? $BaitRegionInfo{$regID}->TDn->tgt : $BaitRegionInfo{$regID}->TDn->max;
        my $regLen = $BaitRegionInfo{$regID}->Len;
        my $spacing = CalculateSpacing($regLen,$density);
        my $numBaits = CalculateNumBaits($regLen,$density,$regID);
        my $spare = $regLen - (($numBaits - 1) * $spacing + $baitLen);
        my $start = $BaitRegionInfo{$regID}->Idx + floor($spare / 2);
        for(my $i = 0; $i < $numBaits; $i++){
            my $header = $regID."+".($start + $i * $spacing);
            my $seq = substr($BaitRegionInfo{$regID}->Seq,$start + $i * $spacing,$baitLen);
            my $seqObj = Bio::Seq->new(-id => $header, seq => $seq);
            print $OUT $seqObj;
        }
    }
}


sub DetermineTilingClass(){
    print STDERR "[DEBUG] DetermineTilingClass\n" if($bDebug);
    my $nOrg = scalar(keys %RegionOrder);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach((values %RegionOrder));
    my $bDone = 0;
    my $MedClassBaits = $maxBaits;
    my @OrgList = sort keys %RegionOrder;
    my @Class = ('Med') x $nOrg;
    my $nMedClass = $nOrg;
    my %BaitCountPerOrg = GetBaitCountsPerOrg();
    #while(my ($org, $countObj) = each %BaitCountPerOrg){
    #    print STDERR "$org: min->$countObj->{min} - max->$countObj->{max}\n";
    #}
    my %SeenIdeal;
    my $Ideal;
    while(!$bDone){
        $bDone = 1;
        if(!$nMedClass){
            warn "[WARNING] No Med Class Organisms: All orgs at either minimum or maximum tiling\n";
            last;
        }
        $Ideal = $MedClassBaits / $nMedClass;
        @Class = ('Med') x $nOrg;
        $nMedClass = $nOrg;
        $MedClassBaits = $maxBaits;
        if ($Ideal <= 0){
            TrimBaits('Min',-$Ideal*$nMedClass+1);
            %BaitCountPerOrg = GetBaitCountsPerOrg();
            $bDone = 0;
        }
        #Use full bait counts assign class to organism,
        #adjust medclassbaits by region. in case of conflict
        #between max and min classes, max wins (max =
        #underbaited org)
        my %MinRegions;
        my %MaxRegions;
        for(my $i = 0; $i < $nOrg; $i++){
            next unless($Class[$i] eq 'Med');
            my $org = $OrgList[$i];
            if ($BaitCountPerOrg{$org}->{max} < $Ideal){
                $Class[$i] = 'Max';
                foreach my $regID (@{$RegionOrder{$org}}){
                    $MaxRegions{$regID} = 1;
                    delete $MinRegions{$regID} if (exists $MinRegions{$regID});
                }
                $nMedClass--;
                $bDone = 0;
            } elsif ($BaitCountPerOrg{$org}->{min} > $Ideal){
                $Class[$i] = 'Min';
                foreach my $regID (@{$RegionOrder{$org}}){
                    $MinRegions{$regID} = 1 unless (exists $MaxRegions{$regID});
                }
                $nMedClass--;
                $bDone = 0;
            }
        }
        foreach my $regID ((keys %MinRegions, keys %MaxRegions)){
            my $density = (exists $MinRegions{$regID}) ? $BaitRegionInfo{$regID}->TDn->min : $BaitRegionInfo{$regID}->TDn->max;
            $MedClassBaits -= CalculateNumBaits($BaitRegionInfo{$regID}->Len,$density,$regID);
        }
        if(exists $SeenIdeal{$Ideal}){
            $bDone = 1;
        } else{
            $SeenIdeal{$Ideal} = 1;
        }
    }
    my %OrgTilingClass;
    @OrgTilingClass{@OrgList} = @Class;
    return ($Ideal,%OrgTilingClass);
}



sub SetRegionTargetTiling(%){
    print STDERR "[DEBUG] SetRegionTargetTiling\n" if($bDebug);
    my %OrgTilingClass = @_;
    my $nOrg = scalar(keys %OrgTilingClass);
    my $MedClassBaits = $maxBaits;
    my $nMedClass=0;
    #Set Class of Regions based on Class of organisms the regions cover
    foreach my $org (sort {scalar(@{$RegionOrder{$a}}) <=> scalar(@{$RegionOrder{$b}})} keys %OrgTilingClass){
        my $class = $OrgTilingClass{$org};
        $nMedClass++ if($class eq 'Med');
        foreach my $regID (@{$RegionOrder{$org}}){
            #Min has lowest precedence
            if(defined $BaitRegionInfo{$regID}->Cls($class)){
                next if($BaitRegionInfo{$regID}->Cls($class) eq 'Max' or $class eq 'Min');
            }
            $BaitRegionInfo{$regID}->Cls($class);
        }
    }
    #Set Tiling of Min and Max Class Regions
    unless($nMedClass == $nOrg){ #All Orgs are Med Class
        foreach my $regID (sort keys %BaitRegionInfo){
            next unless(defined $BaitRegionInfo{$regID}->Cls);
            next if($BaitRegionInfo{$regID}->Cls eq 'Med');
            my $density = $BaitRegionInfo{$regID}->TDn->min;
            $density = $BaitRegionInfo{$regID}->TDn->max if($BaitRegionInfo{$regID}->Cls eq 'Max'); 
            $BaitRegionInfo{$regID}->TDn->tgt($density);
            $MedClassBaits -= CalculateNumBaits($BaitRegionInfo{$regID}->Len,$density,$regID);
        }
    }

    #Set Tiling of Med Class Regions
    if($nMedClass){ #Otherwise all are either min or max and have been completely set
        my $Ideal = $MedClassBaits / $nMedClass;
        foreach my $org (sort {scalar(@{$RegionOrder{$a}}) <=> scalar(@{$RegionOrder{$b}})} keys %OrgTilingClass){
            my $class = $OrgTilingClass{$org};
            next unless($class eq 'Med');
            #ensure Regions are already sorted with highly shared regions last, This should be
            #unnecessary
            for(my $i = 0; $i < @{$RegionOrder{$org}} - 1; $i++){
                if(@{$BaitRegionInfo{$RegionOrder{$org}->[$i]}->Org} > @{$BaitRegionInfo{$RegionOrder{$org}->[$i+1]}->Org}){
                    @{$RegionOrder{$org}} = sort {
                        scalar(@{$BaitRegionInfo{$a}->Org}) <=> scalar(@{$BaitRegionInfo{$b}->Org})
                        } @{$RegionOrder{$org}};
                    last;
                }
            }
            #Sort the region list so that fixed regions come last, but still ensuring that highly
            #shared regions are last, and reverse it so that regions are processed from most to least
            #constrained
            my @regIDList = reverse sort {
                ((defined $BaitRegionInfo{$b}->TDn->tgt) ? 0 : 1) <=> ((defined $BaitRegionInfo{$a}->TDn->tgt) ? 0 : 1)
                } @{$RegionOrder{$org}};
            for(my $i = 0; $i < @regIDList; $i++){

                my $regID = $regIDList[$i];
                next if(defined $BaitRegionInfo{$regID}->TDn->tgt);#Region's density is already fixed
                my $Weight = 0;
                #Set the weight at the maximum across organisms (Better to overbait high, than
                #underbait low
                foreach my $curOrg (@{$BaitRegionInfo{$regID}->Org}){
                    my $curWeight = CalculateWeightofRegionWithinOrg($Ideal,$curOrg,$regID);
                    $Weight = $curWeight if($curWeight > $Weight);
                }
                my $regLen = $BaitRegionInfo{$regID}->Len;
                my $nBaits = (1-$Weight) * CalculateNumBaits($regLen,$BaitRegionInfo{$regID}->TDn->min,$regID);
                $nBaits += $Weight * CalculateNumBaits($regLen,$BaitRegionInfo{$regID}->TDn->max,$regID);
                my $density = CalculateTilingDensity($regLen,$nBaits);
                #print STDERR "$regID\t$Weight\t$nBaits\t$density\n";
                $BaitRegionInfo{$regID}->TDn->tgt(ceil($density));
            }
        }
    }
}

sub CalculateWeightofRegionWithinOrg($$$){
    print STDERR "[DEBUG] CalculateWeightofRegionWithinOrg\n" if($bDebug);
    my($Ideal,$org,$regID) = @_;
    my $distFromIdeal = $Ideal;
    my $potential = 0;
    foreach my $otherRegID (@{$RegionOrder{$org}}){
        #next if($otherRegID eq $regID);
        if(defined $BaitRegionInfo{$otherRegID}->TDn->tgt){
            my $density = $BaitRegionInfo{$otherRegID}->TDn->tgt;
            $distFromIdeal -= CalculateNumBaits($BaitRegionInfo{$otherRegID}->Len,$density,$regID);
        } else {
            my $density = $BaitRegionInfo{$otherRegID}->TDn->max;
                $density += $BaitRegionInfo{$otherRegID}->TDn->min;
                $density = floor($density/2);
                $density = $BaitRegionInfo{$otherRegID}->TDn->max;
            my $len = $BaitRegionInfo{$otherRegID}->Len;
            $potential += CalculateNumBaits($len,$density,$regID); 
        }
    }
    my $Weight = 0;
    if($distFromIdeal > 0){
        $Weight = 1 - ($potential / ($distFromIdeal + $potential));
    }
    return $Weight;
}


sub CalculateTilingDensity($$){
    print STDERR "[DEBUG] CalculateTilingDensity\n" if($bDebug);
    my ($regLen,$nBaits) = @_;
    if($nBaits + 1 == ceil(($regLen+1)/$baitLen)){
        return 1;
    } elsif($regLen < 2*$baitLen){
        return floor($nBaits*($regLen - $baitLen)/($regLen - $baitLen + 1))+1;
    } else {
        return floor($nBaits*($baitLen - 1)/($regLen - $baitLen + 1))+1;
    }
}

sub ApplyMaxBaitPerOrgFilter(;$){
    print STDERR "[DEBUG] ApplyMaxBaitPerOrgFilter\n" if($bDebug);
    my $bQuiet = shift(@_);
    $bQuiet = 0 unless($bQuiet);
    my %BaitCountsPerOrg = GetBaitCountsPerOrg();
    my %nTrimPerOrg;
    my $sumCount = 0;
    foreach my $org (keys %BaitCountsPerOrg){
        $nTrimPerOrg{$org}=0;
        my $count = $BaitCountsPerOrg{$org}->{max};
        $sumCount += $count;
        $nTrimPerOrg{$org} = $count - $maxBaitsPerOrg if($count > $maxBaitsPerOrg);
    }
    TrimBaitsFromOrganisms('Max',($bQuiet) ? $sumCount: undef,%nTrimPerOrg);
}

sub TrimBaits($$){
    print STDERR "[DEBUG] TrimBaits\n" if($bDebug);
    my ($densMode,$nBaitsToTrim) = (@_);
    my $TotalTrimmed = 0;
    my $bCycle = 0;
    my $remaining = $nBaitsToTrim;
    while(!$bCycle and $remaining > 0){ #Repeat until the desired number of baits have been removed or no more baits can be removed
        $bCycle = 1;
        my %BaitCountsPerOrg = GetBaitCountsPerOrg();
        $BaitCountsPerOrg{$_} = $BaitCountsPerOrg{$_}->{lc($densMode)} foreach(keys %BaitCountsPerOrg);
        my @Orgs = sort {$BaitCountsPerOrg{$b} <=> $BaitCountsPerOrg{$a}} (keys %BaitCountsPerOrg);
        my %nTrimPerOrg;
        @nTrimPerOrg{@Orgs} = (0) x @Orgs;
        for(my $i = 0; $i < @Orgs; $i++){
            my $modulus = 0;
            my $diff = 0;
            if($i < $#Orgs and $BaitCountsPerOrg{$Orgs[$i+1]} > $minBaitsPerOrg){
                $diff = $BaitCountsPerOrg{$Orgs[$i]} - $BaitCountsPerOrg{$Orgs[$i + 1]};
            } else {
                $diff = $BaitCountsPerOrg{$Orgs[$i]} - $minBaitsPerOrg;
            }
            $diff = 0 if($diff < 0);
            if($diff * ($i + 1) > $remaining){
                $modulus = $remaining % ($i+1);
                $diff = int($remaining/($i+1));
            }
            for(my $j = 0; $j <= $i; $j++){
                my $maxTrim = $BaitCountsPerOrg{$Orgs[$j]} - $minBaitsPerOrg;
                $maxTrim = ($maxTrim > 0) ? $maxTrim : 0;
                my $availableTrim = ($maxTrim > $nTrimPerOrg{$Orgs[$j]}) ? $maxTrim - $nTrimPerOrg{$Orgs[$j]} :0;
                my $ToTrim = ($availableTrim > $diff + $modulus) ? ($diff + $modulus) : $availableTrim;
                #print STDERR "[$i,$j] $availableTrim($nTrimPerOrg{$Orgs[$j]},$maxTrim)\t$diff($BaitCountsPerOrg{$Orgs[$i]} - $BaitCountsPerOrg{$Orgs[$i+1]})+$modulus\t$remaining\t=$ToTrim\n";
                $modulus = ($diff + $modulus) - $ToTrim;
                $nTrimPerOrg{$Orgs[$j]} += $ToTrim;
                $remaining -= $ToTrim;
                last if(!$remaining);
            }
            if(!$remaining){
                $bCycle = 0;
                last;
            }
        }
        my %Val = (max => -1, min => 99**99, sum => 0);
        foreach my $org (keys %nTrimPerOrg){
            $Val{sum} += $nTrimPerOrg{$org};
            $Val{max} = $nTrimPerOrg{$org} if($Val{max} < $nTrimPerOrg{$org});
            $Val{min} = $nTrimPerOrg{$org} if($Val{min} > $nTrimPerOrg{$org});
        }
        my $trimmed = TrimBaitsFromOrganisms($densMode,$nBaitsToTrim,%nTrimPerOrg);
        $TotalTrimmed += $trimmed;
        $remaining = $nBaitsToTrim - $TotalTrimmed;
        $bCycle = 1 if($trimmed == 0);
        #print STDERR "$Val{sum} ($Val{min} -> $Val{max})\t$remaining\t$TotalTrimmed\n";
    }
    if($TotalTrimmed < $nBaitsToTrim){
        warn "[WARNING] Could not trim $nBaitsToTrim baits due to region sharing or size: Trimmed $TotalTrimmed baits\n";
    }
}

sub TrimBaitsFromOrganisms($$%){
    #Entire Regions may be trimmed away
    #As The first 3 sorting priorities are unaffected by region length, and ties past that point are
    #unlikely, Baits are successively removed from the lowest priority regions
    print STDERR "[DEBUG] TrimBaitsFromOrganism\n" if($bDebug);
    my $densMode = shift(@_);
    my $trimCap = shift(@_);
    my %nBaitsPerOrg = @_;
    my %Record = @_;
    my $bTrimCap = 1;
    unless(defined $trimCap){
        $bTrimCap = 0;
        $trimCap = 0;
        $trimCap += $nBaitsPerOrg{$_} foreach(keys %nBaitsPerOrg);
    }
    my %UntrimmableRegions;
    my $TotalTrimmed = 0;
    foreach my $org (sort {$nBaitsPerOrg{$b} <=> $nBaitsPerOrg{$a}} keys %nBaitsPerOrg){
        next unless($nBaitsPerOrg{$org});
        my $remaining = $nBaitsPerOrg{$org};
        foreach my $regID (reverse @{$RegionOrder{$org}}){
            next if(exists $UntrimmableRegions{$regID});
            my $toTrim = 0;
            my $regLen = $BaitRegionInfo{$regID}->Len;
            my $density = GetRegionDensityByMode($densMode,$regID);
            my $availableToTrim = CalculateNumBaits($regLen,$density,$regID);
            my $init = $availableToTrim;
            my $maxTrim = $availableToTrim;
            foreach my $oOrg (@{$BaitRegionInfo{$regID}->Org}){
                if($nBaitsPerOrg{$oOrg} < $availableToTrim){
                    $availableToTrim = $nBaitsPerOrg{$oOrg};
                }
            }
            if(!$availableToTrim){ #This region cannot and never will be trimmed
                $UntrimmableRegions{$regID} = 1;
                next;
            }
            my $bDone = 0;
            my $Trimmed = ($availableToTrim < $remaining) ? $availableToTrim : $remaining;
            $Trimmed = $trimCap - $TotalTrimmed if($TotalTrimmed + $Trimmed > $trimCap);
            #print STDERR "$regID: $Trimmed(av:$availableToTrim(dn:$density, rL: $regLen), rm:$remaining, tc: $trimCap)\n";
            if($Trimmed == $maxTrim){
                RemoveRegion($regID);
            } else {
                TrimBaitsFromRegion($regID,$Trimmed,$density);
            }
            $remaining -= $Trimmed;
            $TotalTrimmed += $Trimmed;
            foreach my $oOrg (@{$BaitRegionInfo{$regID}->Org}){
                $nBaitsPerOrg{$oOrg} -= $Trimmed;
            }
            last if(!$remaining);
            last if($TotalTrimmed == $trimCap)
        }
        if($nBaitsPerOrg{$org} && !$bTrimCap){
            warn "[WARNING] Could not trim all $Record{$org} baits from $org due to excessive sharing:".
                " Trimmed @{[$Record{$org} - $nBaitsPerOrg{$org}]}\n";
        }
    }
    return $TotalTrimmed;
}

sub RemoveRegion($){
    print STDERR "[DEBUG] RemoveRegion\n" if($bDebug);
    my ($regID) = @_;
    #To Remove A Region, all references to in in the Region Order structure need to be removed,
    #As eveything is done in reference to this structure, it needn't be removed from the
    #BaitRegionInfo Structure
    foreach my $org (@{$BaitRegionInfo{$regID}->Org}){
        my $i = $#{$RegionOrder{$org}};
        $i-- until($RegionOrder{$org}->[$i] eq $regID);
        splice(@{$RegionOrder{$org}},$i,1);
    }
}

sub TrimBaitsFromRegion($$$){
    print STDERR "[DEBUG] TrimBaitsFromRegion\n" if($bDebug);
    my ($regID,$nBaits,$density) = @_;
    #Trimming a Region invloves reducing its nominal length,
    #   adjusting its start position so it remains centered,
    #   and updating its tiling object
    my $bDone = 0;
    my $curBaits = CalculateNumBaits($BaitRegionInfo{$regID}->Len,$density,$regID);
    my $targetBaits = $curBaits - $nBaits;
    if($targetBaits < 1){
        RemoveRegion($regID);
        return;
    }
    while(!$bDone){
        $bDone = 1;
        my $regLen = $BaitRegionInfo{$regID}->Len;
        my $targetSpacing = CalculateSpacing($regLen,$density);
        my $newLen = $targetBaits * $targetSpacing + $baitLen -1;
        $BaitRegionInfo{$regID}->Idx(floor(($regLen - $newLen)/2));
        $BaitRegionInfo{$regID}->Len($newLen);
        my $Tdn = GetTilingDens($newLen);
        if(defined $BaitRegionInfo{$regID}->TDn->tgt){
            $Tdn->tgt($BaitRegionInfo{$regID}->TDn->tgt);
            $Tdn->tgt($Tdn->max) if($Tdn->tgt > $Tdn->max);
        }
        $density = $Tdn->max if($density > $Tdn->max);
        $curBaits = CalculateNumBaits($newLen,$density,$regID);
        $BaitRegionInfo{$regID}->TDn($Tdn);
        $bDone = 0 if($curBaits > $targetBaits);
    }
}


sub OptimizeMaxBaitsPerOrg(){
    print STDERR "[DEBUG] OptimizeMaxBaitsPerOrg\n" if($bDebug);
    my $nMinClass = 2;
    my @MinOrgs;
    for(my $i = 0; $i < $nMinClass; $i++){
        $nMinClass = 0;
        my (undef,%OrgTilingClass) = DetermineTilingClass();
        my %BaitCountPerOrg = GetBaitCountsPerOrg();
        $BaitCountPerOrg{$_} = $BaitCountPerOrg{$_}->{max} foreach(keys %BaitCountPerOrg);
        my $MaxMid = 0;

        my %MinOrgSet;
        @MinOrgSet{@MinOrgs} = (1) x @MinOrgs;
        while(my ($org,$class) =each %OrgTilingClass){
            $MaxMid = $BaitCountPerOrg{$org} if($class eq 'Med' and $BaitCountPerOrg{$org} > $MaxMid);
            next unless($class eq 'Min');
            $nMinClass++;
            push(@MinOrgs,$org) unless(exists $MinOrgSet{$org});
        }
        if($nMinClass > 0){
            if($i == 0){
                @MinOrgs = sort {$BaitCountPerOrg{$b} <=> $BaitCountPerOrg{$a}} @MinOrgs;
            }
            if($i >= $nMinClass - 2){
                $maxBaitsPerOrg = $MaxMid;#$BaitCountPerOrg{$MinOrgs[$i]};
            } else {
                $maxBaitsPerOrg = $BaitCountPerOrg{$MinOrgs[$i+1]};
            }
            ApplyMaxBaitPerOrgFilter(1);
        } elsif(!$maxBaitsPerOrg){
            $maxBaitsPerOrg = $maxBaits;
        }
        PrintUpdate("Optimization Iteration @{[$i+1]}: Max @ $maxBaitsPerOrg",1) if($bDebug);
    }
}

sub OptimizeMaxRegionsPerOrg(){
    print STDERR "[DEBUG] OptimizeMaxRegionsPerOrg\n" if($bDebug);
    my $nMinClass = 3;
    my @MinOrgs;
    for(my $i = 0; $i < $nMinClass; $i++){
        $nMinClass = 0;
        my (undef,%OrgTilingClass) = DetermineTilingClass();
        my $MaxMid = 0;
        my %MinOrgSet;
        @MinOrgSet{@MinOrgs} = (1) x @MinOrgs;
        while(my ($org,$class) =each %OrgTilingClass){
            if($class eq 'Med' and scalar(@{$RegionOrder{$org}}) > $MaxMid){
                $MaxMid = scalar(@{$RegionOrder{$org}});
            }
            next unless($class eq 'Min');
            $nMinClass++;
            push(@MinOrgs,$org) unless(exists $MinOrgSet{$org});
        }
        if($nMinClass > 0){
            if($i == 0){
                @MinOrgs = sort {scalar(@{$RegionOrder{$b}}) <=> scalar(@{$RegionOrder{$a}})} @MinOrgs;
            }
            if($i == $nMinClass - 1){
                $maxRegPerOrg = $MaxMid;
            } else {
                #print STDERR $nMinClass,"_",scalar(@MinOrgs),"&",$i,"*",$MinOrgs[$i],"\n";
                $maxRegPerOrg = scalar(@{$RegionOrder{$MinOrgs[$i]}});
            }
            ApplyMaxRegionFilter(%RegionOrder,1);
        } elsif(!$maxRegPerOrg){
            $maxRegPerOrg = $maxBaits;
        }
        PrintUpdate("Optimization Iteration @{[$i+1]}: Max @ $maxRegPerOrg",1) if($bDebug);
    }
}

sub PrintUpdate($$){
    my ($message,$bTime) = @_;
    my $string = "[INFO] $message";
    my $stepTime = time;
    $string .= " (".($stepTime - $LastTime)."s)" if($bTime);
    print STDERR $string,"\n";
    $LastTime = $stepTime;
}
