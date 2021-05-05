#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::TreeIO;
use Scalar::Util qw(looks_like_number);
use Getopt::Std;
use Class::Struct;
use POSIX;
use sort 'stable';
use FindBin;
use File::Basename;
use lib File::Spec->catdir($FindBin::RealBin,'..','lib');
use HUBDesign::Util qw(ProcessNumericOption OpenFileHandle LoadConfig);
use HUBDesign::Logger;


#============================================================
#Declarations


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
#tid = The actual id of the taxon target be this region
#Ex. A 400bp Sequence which has been trimmed to 4 baits at 4x coverage would have: Len=170, Idx=115
#Ex. A 150bp Sequence which has been trimmed to 37 baits at 35x coverage would have Len=149, Idx=0
struct(REGION => {Org => '@', Tgt => '$', Seq => '$', Pos => '$',
        TDn => 'TILING', Cls => '$', Idx => '$', Len => '$',tid => '$'});
sub PrintREGION($$); #Usage RegObj->PrintREGION($file Handle);


sub ProcessDelimStr_ListFile($$$$);
sub LoadTree($);	                    #Usage my $TreeObj = LoadTree($file);
sub LoadRegionInfo($$);                     #Usage my %RegionInfo = LoadRegionInfo($file,$TreeObj);
sub GetDescendentLeaves($$);	            #Usage my @IdList = GetDescendentLeaves($taxon_id,$tree)
sub GetSortedRegions();
sub sortMaxSpread(@);
sub MinDistance($@);
sub GetTilingDens($);                       #Usage my TilingObj = GetTiling(seqLen);
sub CalculateSpacing($$);
sub CalculateNumBaits($$;$);
sub GetBaitCount($@);                       #Usage sub GetMaximumBaitCount(('Min','Max','Tgt',$Int),@ListOfOrderedRegIDListRefs);
sub GetRegionDensityByMode($$);             #Usage my $density = GetRegionDensityByMode($Mode,$regID);
sub ApplyMaxRegionFilter(\%;$);
sub OutputBaits(@);                         #Usage OutputBaits(@ListOfOrderedRegIDListRefs)
sub DetermineTilingClass();                 #Usage my @OrgTilingClass = DetermineTilingClass();
sub GetBaitCountsPerOrg();
sub SetRegionTargetTiling(%);               #Usage SetRegionTargetTiling(%OrgTilingClass);
sub CalculateWeightofRegionWithinOrg($$$);  #Usage fncname($IdealBaitsPerOrg,$orgID,$regID);
sub CalculateTilingDensity($$);             #Usage fncname($regLen,$numBaits);
sub TrimBaits($$);                          #Usage $BaitsTrimmed = TrimBaits($NumBaitstoTrim);
sub TrimBaitsFromOrganisms($$%);            #Usage $BaitsTrimmed = TrimBaitsFromOrganism($densityMode,(Org1 => NumToRemove))
sub ApplyMaxBaitPerOrgFilter(;$);           
sub OptimizeMaxBaitsPerOrg();               
sub OptimizeMaxRegionsPerOrg();
sub ParseConfig($);	                    #Usage ParseConfig($configFile);

my $MIN_BAITS_PER_ORG = 30;
my $_ID_PADDING = 6;
my %DEFAULT = (b => 0, d => 1, i => "BaitInfo.tsv", l => 75, n => 0, r => 0);


#============================================================
#Handling Input

my %opts;
getopts('DOb:d:i:l:n:p:r:C:aPSvh',\%opts);
ParseConfig($opts{C});

if(exists $opts{h} or @ARGV < 2 or !exists $opts{n}){
    my $Usage = "@{[basename($0)]} [-options] -n maxBaits Tree RegionInfo > TiledBaits.fna";
    if(exists $opts{h}){
        die "===Description\n".
            "\tGiven a set of Bait Regions and some metadata, outputs a set of baits which\n".
            "\t\toptimally tile the regions.\n".
            "\tOptimal is baits centered within regions, with maximally balanced baits per species\n".
            "\t\twith maximally spaced small regions given priority\n".
            "===Usage\n".
            "\t$Usage\n".
            "===Required Inputs\n".
            "\tTree\tThe tree detailing the relationships between taxon ids\n".
            "\tRegionInfo\tA tab delimited file with 5 columns:\n".
            "\t\tRegion ID, Taxon ID, Cluster ID, Position within Cluster, and Sequence\n".
           "\t-n INT[0,∞)\tThe maximum number of baits in the output[Default: $DEFAULT{n}]\n".
           "===Options\n".
           "\t-b INT[0,∞)\tThe maximum number of baits for any target organism [Default: $DEFAULT{b}]\n".
           "\t-d INT[1,l]\tThe minimum tiling density for any bait region [Default: $DEFAULT{d}x]\n". 
           "\t-i PATH\tThe file to which the bait metadata will be saved (Default: $DEFAULT{i})\n".
           "\t-l INT[1,∞)\tThe length of baits to output [Default: $DEFAULT{l}bp]\n".
           "\t-p string\tA comma separated string or list file detailing the order of prescedence\n".
           "\t\tfor clusters. By default all are considered equal\n".
           "\t-r INT[0,∞)\tThe maximum number of regions for any target organism [Default: $DEFAULT{r}]\n".
           "\t-C PATH\t Path to a HUBDesign Config file (Default: $FindBin::RealBin/../HUBDesign.cfg)\n".
           "\tNote: For b, n, and r a value of zero indicates no maximum\n".
           "===Flags\n".
           "\t-a\tAuto-detect max regions and baits per Organism to improve Proportional Tiling\n".
           "\t-P\tDisable Proportional tiling, and instead find a maximum constant tiling\n".
           "\t-S\tSkip Region Sorting\n".
           "\t-v\tVerbose Output\n".
           "\t-h\tDisplay this message and exit\n";
    } else {
        die "Usage: $Usage\n\tUse -h for more info\n";
    }
}

foreach (keys %DEFAULT){
    $opts{$_} = exists $opts{$_} ? $opts{$_} : $DEFAULT{$_};
}

#Init Logger Verbosity
$opts{v} = (exists $opts{v}) ? "INFO" : "WARNING";
$opts{v} = "DEBUG" if(exists $opts{D});
my $Logger = HUBDesign::Logger->new(level => $opts{v});


my %FileDict;
@FileDict{qw(Tree BRInfo BaitInfo)} = (@ARGV[0 .. 1], $opts{i});
#Check if given files exist
while (my ($type,$file) = each %FileDict){
    next if($type eq "BaitInfo"); 
    $Logger->Log("$type file ($file}) doesn't exist","ERROR") unless(-e $file);
}



## Process Simple Numeric Options
my $maxBaitsPerOrg = ProcessNumericOption($opts{b},0,0,undef,1,"Maximum Baits per Organism");
my $baitLen = ProcessNumericOption($opts{l},75,1,undef,1,"BaitLength");
my $minDensity = ProcessNumericOption($opts{d},1,1,$baitLen,1,"Minimum Tiling Density");
my $maxBaits = ProcessNumericOption($opts{n},0,0,undef,1,"Maximum Baits");
my $maxRegPerOrg = ProcessNumericOption($opts{r},0,0,undef,1,"Maximum Regions per Organism");
## Handle Cluster Order
my @TargetSeqOrder = ProcessDelimStr_ListFile($opts{p},",",[],"Target Sequence Order");
if(exists $opts{a}){
    if(exists $opts{P}){
        $Logger->Log("Auto-detect option is useless if proportional tiling is disabled","WARNING");
    }
    if($maxBaitsPerOrg){
        $Logger->Log("Auto-detect overuling specified maxBaitsPerOrg","WARNING");
        $maxBaitsPerOrg = 0;
    }
    if($maxRegPerOrg){
        $Logger->Log("Auto-detect overuling specified maxRegPerOrg","WARNING");
        $maxRegPerOrg = 0;
    }
}
#Handle Flags
my $bRegSort = exists $opts{S} ? 0 : 1;
my $bDebug = (exists $opts{D}) ? 1 : 0;
my $OutMode = ('TaxaInfo','Baits')[exists $opts{O} ? 0 : 1];

#Log Parameters
my %params = (  maxBaitsPerOrg => ($maxBaitsPerOrg) ? $maxBaitsPerOrg : "∞",
                baitLen => $baitLen, minDensity => $minDensity,
                maxBaits => $maxBaits, 
                maxRegPerOrg => ($maxRegPerOrg) ? $maxRegPerOrg : "∞",
                TargetSeqOrder => (@TargetSeqOrder) ? join(",",@TargetSeqOrder) : "No Order",
                TilingMode => (exists $opts{P}) ? 'Maximal' : 'Proportional',
                RegionInfoFile => $FileDict{BRInfo},
                BaitInfoFile => $FileDict{BaitInfo},
                TreeFile => $FileDict{Tree});
@params{qw(maxBaitsPerOrg maxRegPerOrg)} = ("Auto-Detect") x 2 if(exists $opts{a});
$Logger->LogParameters(map {($_,$params{$_})} (sort keys %params));

#============================================================
#Main

my $TreeObj = LoadTree($FileDict{Tree});
my %BaitRegionInfo = LoadRegionInfo($FileDict{BRInfo},$TreeObj);
my %RegionOrder = GetSortedRegions();
ApplyMaxRegionFilter(%RegionOrder) if($maxRegPerOrg);
my $AbsMaxBaits = GetBaitCount('Max',(values %RegionOrder));
ApplyMaxBaitPerOrgFilter() if($maxBaitsPerOrg);

my %OrgTilingClass;
if($maxBaits && $AbsMaxBaits > $maxBaits){
    if(exists $opts{P}){ #Select maximal tiling keeping # Baits below thresh
        $Logger->Log("Determining Maximal Constant Tiling Density","INFO");
        my $density = $minDensity;
        $density++ while(GetBaitCount($density+1,(values %RegionOrder)) <= $maxBaits);
        foreach my $regID (sort keys %BaitRegionInfo){
            my $min = $BaitRegionInfo{$regID}->TDn->min;
            my $max = $BaitRegionInfo{$regID}->TDn->max;
            my $tgt = ($density < $min) ? $min : $density;
            $tgt = ($tgt > $max) ? $max : $tgt;
            $BaitRegionInfo{$regID}->TDn->tgt($tgt);
        }
        $Logger->Log("Found Maximal Constant Tiling Density at ${density}x","INFO");
        (undef,%OrgTilingClass) = DetermineTilingClass();
    } else { #Scale Tiling to balance #Baits per Org
         if(exists $opts{a}){
             OptimizeMaxRegionsPerOrg();
             OptimizeMaxBaitsPerOrg();
         }
         $Logger->Log("Selecting proportional tiling densities...","INFO");
         my $Ideal;
         ($Ideal,%OrgTilingClass) = DetermineTilingClass();
         SetRegionTargetTiling(%OrgTilingClass);
         if(exists $opts{v}){
             $Ideal = sprintf("%.0f",$Ideal);
             $Logger->Log("Balanced tiling densities target $Ideal baits per organism","INFO");
         }
    }
    my $curBaits = GetBaitCount('Tgt',(values %RegionOrder));
    while($curBaits > $maxBaits){
        $Logger->Log("Currently have $curBaits: Trimming baits down to $maxBaits...","INFO");
        TrimBaits('Tgt',$curBaits - $maxBaits);
        my $nowBaits = GetBaitCount('Tgt',(values %RegionOrder));
        if(exists $opts{v}){
            $Logger->Log("Trimmed @{[$curBaits - $nowBaits]} baits","INFO");
        }
        $curBaits = $nowBaits;
    }
} elsif($maxBaits) {
    $Logger->Log("Maximum Possible Baits is ≤ specified maximum : All baits to be output at maximum density","WARNING");
    foreach my $RegObj (values %BaitRegionInfo){
        $RegObj->TDn->tgt($RegObj->TDn->max);
    }
}
if($OutMode eq 'Baits'){
    OutputBaits((values %RegionOrder));
} elsif($OutMode eq 'TaxaInfo') {
    my %BaitCountPerOrg = GetBaitCountsPerOrg();
    foreach my $org (sort {$BaitCountPerOrg{$a}->{tgt} <=> $BaitCountPerOrg{$b}->{tgt}} keys %BaitCountPerOrg){
        print "$org\t$BaitCountPerOrg{$org}->{tgt}\t$OrgTilingClass{$org}\n";
    }
}

$Logger->Log("Done","INFO");

#============================================================
#Subroutines

sub ProcessDelimStr_ListFile($$$$){
    $Logger->Log("ProcessDelimStr_ListFile","DEBUG");
    my($val,$delim,$default,$varName) = @_;
    return @{$default} unless(defined $val);
    my @valList = split($delim,$val);
    my @Output;
    foreach my $file (@valList){
        if(-e $file){
            my $fh = OpenFileHandle($file,$varName,"ERROR");
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
    $Logger->Log("PrintREGION","DEBUG");
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


sub LoadTree($){
    my $file = shift;
    my $TreeIOObj = Bio::TreeIO->new(-file => $file);
    return $TreeIOObj->next_tree;
}

sub LoadRegionInfo($$){
    my ($file,$treeObj) = @_;
    $Logger->Log("Loading bait regions...","INFO");
    my %BRInfo; # Hash with keys of region ids and values of REGION
    my $fh = OpenFileHandle($file,"RegionInfo");
    while(my $line = <$fh>){
        chomp($line);
        my ($id,$taxon_id,$target,$pos,$seq) = split(/\t/,$line);
        if(exists $BRInfo{$id}){
            $Logger->Log("Ignoring additional entries for region $id on line $.","WARNING");
            next;
        }
        my @orgList = GetDescendentLeaves($taxon_id,$treeObj);
        @orgList = ($taxon_id) unless (@orgList);
        $BRInfo{$id} = REGION->new(Org => \@orgList,Tgt => $target,Pos=>$pos,Idx=>0,Seq=>$seq,
            Len => length($seq),TDn => GetTilingDens(length($seq)), tid => $taxon_id);
    }
    close($fh);
    $Logger->Log("Loaded @{[scalar(keys %BRInfo)]} bait regions","INFO");
    return %BRInfo;
}

sub GetDescendentLeaves($$){
    my ($taxon_id,$treeObj) = @_;
    my @leaves;
    my ($node) = $treeObj->find_node(-id => $taxon_id);
    return () unless (defined $node);
    foreach my $child ($node->get_all_Descendents){
       push(@leaves,$child->id) if($child->is_Leaf); 
    }
    return @leaves;
}

sub GetSortedRegions() {
    $Logger->Log("Sorting regions within organisms..","INFO");
    my %RegOrder;
    ##Assign each region to its organism
    foreach my $regID (sort keys %BaitRegionInfo){
        my $RegionObj = $BaitRegionInfo{$regID};
        foreach my $org (@{$RegionObj->Org}){
            unless(exists $RegOrder{$org}){
                $RegOrder{$org} = [];
            }
            push(@{$RegOrder{$org}},$regID);
        }
    }
    return %RegOrder unless($bRegSort);
    my %TargetSeqRank;
    if(@TargetSeqOrder){
        for(my $i = 0; $i < @TargetSeqOrder; $i++){
            $TargetSeqRank{$TargetSeqOrder[$i]} = $i;
        }
    }
    my $OrgCounter = 0;
    my $OrgCount = scalar(keys %RegOrder);
    my $LastTime = time;
    foreach my $org (sort keys %RegOrder){
        my $idListRef = $RegOrder{$org};
        my $elapsed = time - $LastTime;
        $OrgCounter++;
        if($elapsed > 10){
            $Logger->Log("Now sorting $org ($OrgCounter/$OrgCount)","INFO");
        }
        ##Sort first by Size - Prioritizing Small Regions
        @{$idListRef} = sort {$BaitRegionInfo{$a}->Len <=> $BaitRegionInfo{$b}->Len} @{$idListRef};
        ##Sort such that maximal spacing is maintained (Stable sort: Size to break ties)
        my @idxPosList;
        for(my $i = 0; $i < @{$idListRef}; $i++){
            push(@idxPosList,{idx => $i, pos => $BaitRegionInfo{$idListRef->[$i]}});
        }
        @idxPosList = sortMaxSpread(@idxPosList);
        @{$idListRef} = (@{$idListRef})[map {$_->{idx}} @idxPosList];
        ###Sort in order of target precedence (Stable sort : Spacing then Size to break ties)
        if(@TargetSeqOrder){
            foreach  my $regID (@{$idListRef}){
                unless(exists $TargetSeqRank{$BaitRegionInfo{$regID}->Tgt}){
                    $TargetSeqRank{$BaitRegionInfo{$regID}->Tgt} = scalar(@TargetSeqOrder) + 1;
                }
            }
            @{$idListRef} = sort {$TargetSeqRank{$BaitRegionInfo{$a}->Tgt} <=> $TargetSeqRank{$BaitRegionInfo{$b}->Tgt}} @{$idListRef}
        }; 
        ##Sort in order of specificiy - Prioritizing fewer target organisms (Stable Sort)
        @{$idListRef} = sort {scalar(@{$BaitRegionInfo{$a}->Org}) <=> scalar(@{$BaitRegionInfo{$b}->Org})} @{$idListRef};
    }
    return %RegOrder;
}

sub sortMaxSpread(@){
    $Logger->Log("sortMaxSpread","DEBUG");
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

sub GetTilingDens($){
    $Logger->Log("GetTilingDens","DEBUG");
    my ($regLen) = @_;
    my $maxDen = ($regLen >= 2*$baitLen-1) ? $baitLen : $regLen-$baitLen+1;
    my $minDen = ($maxDen < $minDensity) ? $maxDen : $minDensity;
    return TILING->new(max => $maxDen, min => $minDen);
}

sub CalculateSpacing($$){
    $Logger->Log("CalculateSpacing(@_) @ @{[ caller ]}","DEBUG");
    my ($regLen,$density) = @_;
    return $baitLen if ($density == 1);
    return ($regLen < 2*$baitLen) ? floor(($regLen - $baitLen)/($density-1)) : floor(($baitLen - 1)/($density - 1));

}

sub CalculateNumBaits($$;$){
    $Logger->Log("CalculateNumBaits(@_) @ @{[ caller ]}","DEBUG");
    my ($regLen,$density,$regID) = @_;
    my $spacing = CalculateSpacing($regLen,$density);
    if(!$spacing){
        $Logger->Log("ID: $regID","DEBUG");
        PrintREGION($BaitRegionInfo{$regID},*STDERR);
    }
    return ceil(($regLen - $baitLen + 1) / $spacing);
}

sub GetBaitCount($@){
    $Logger->Log("GetBaitCount(@_) @ @{[ caller ]}","DEBUG");
    my $densMode = shift(@_);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(@_);
    my $count = 0;
    foreach my $regID (sort keys %UniqRegID){
        my $density = GetRegionDensityByMode($densMode,$regID); 
        $count += CalculateNumBaits($BaitRegionInfo{$regID}->Len,$density,$regID);
    }
    return $count;
}

sub GetBaitCountsPerOrg(){
    $Logger->Log("GetBaitCountsPerOrg","DEBUG");
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
    $Logger->Log("GetRegionDensityByMode(@_) @{[ caller ]}","DEBUG");
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
        $Logger->Log("Unknown density Mode: $densMode","ERROR");
    }
    return $density;
}

sub ApplyMaxRegionFilter(\%;$){
    my($RegionOrderRef,$bQuiet) = @_;
    $bQuiet = 0 unless(defined $bQuiet);
    my $loglevel = ($bQuiet) ? "DEBUG" : "INFO";
    $Logger->Log("Applying Maximum Regions Per Organism Filter",$loglevel);
    #maxRegPerOrg is assumed to be non-zero
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
            my $message = "Could not acheive at most $maxRegPerOrg regions for $org, due to region sharing: acheived @{[$maxRegPerOrg + $OrgAboveThresh{$org}]} regions";
            $Logger->Log($message,"WARNING");
        }
    }
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(sort (values %RegionOrder));
    $Logger->Log("Filtered out @{[scalar(keys %BaitRegionInfo) - scalar(keys %UniqRegID)]} Regions",$loglevel);
}

sub OutputBaits(@){
    my $BaitCount = GetBaitCount('Tgt',(values %RegionOrder));
    $Logger->Log("Outputting $BaitCount Baits","INFO");
    my $infofh = OpenFileHandle(">$FileDict{BaitInfo}","BaitInfo");
    print $infofh "ProbeID\tTaxonID\tClustID\tClustPos\n";
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(@_);
    my $OUT = Bio::SeqIO->newFh(-format => "fasta");
    my $nextid = 0;
    foreach my $regID (sort keys %UniqRegID){
        my $density = (defined $BaitRegionInfo{$regID}->TDn->tgt) ? $BaitRegionInfo{$regID}->TDn->tgt : $BaitRegionInfo{$regID}->TDn->max;
        my $regLen = $BaitRegionInfo{$regID}->Len;
        my $spacing = CalculateSpacing($regLen,$density);
        my $numBaits = CalculateNumBaits($regLen,$density,$regID);
        my $spare = $regLen - (($numBaits - 1) * $spacing + $baitLen);
        my $start = $BaitRegionInfo{$regID}->Idx + floor($spare / 2);
        for(my $i = 0; $i < $numBaits; $i++){
            my $offset = $start + $i * $spacing;
            my $probe_id = sprintf("Probe_%0*d",$_ID_PADDING,$nextid++);
            my $seq = substr($BaitRegionInfo{$regID}->Seq,$offset,$baitLen);
            my $seqObj = Bio::Seq->new(-id => $probe_id, seq => $seq);
            print $OUT $seqObj;
            print $infofh join("\t",($probe_id,$BaitRegionInfo{$regID}->tid,$BaitRegionInfo{$regID}->Tgt,$BaitRegionInfo{$regID}->Pos + $offset)),"\n";
        }
    }
    close($infofh);
}


sub DetermineTilingClass(){
    $Logger->Log("DetermineTilingClass","DEBUG");
    my $nOrg = scalar(keys %RegionOrder);
    my %UniqRegID;
    @UniqRegID{@{$_}} = (1) x @{$_} foreach(sort (values %RegionOrder));
    my $bDone = 0;
    my $MedClassBaits = $maxBaits;
    my @OrgList = sort keys %RegionOrder;
    my @Class = ('Med') x $nOrg;
    my $nMedClass = $nOrg;
    my %BaitCountPerOrg = GetBaitCountsPerOrg();
    my %SeenIdeal;
    my $Ideal;
    while(!$bDone){
        $bDone = 1;
        if(!$nMedClass){
            $Logger->Log("No Med Class Organisms: All orgs at either minimum or maximum tiling; maxBaits may be too low, or minDensity too high","WARNING");
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
        foreach my $regID ((sort keys %MinRegions, sort keys %MaxRegions)){
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
    $Logger->Log("SetRegionTargetTiling","DEBUG");
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
                $BaitRegionInfo{$regID}->TDn->tgt(ceil($density));
            }
        }
    }
}

sub CalculateWeightofRegionWithinOrg($$$){
    $Logger->Log("CalculateWeightofRegionWithinOrg","DEBUG");
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
    $Logger->Log("CalculateTilingDensity","DEBUG");
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
    my $bQuiet = shift(@_);
    $bQuiet = 0 unless($bQuiet);
    my $logLevel = ($bQuiet) ? "DEBUG" : "INFO";
    $Logger->Log("Applying Maximum Bait Per Organism Filter",$logLevel);
    my %BaitCountsPerOrg = GetBaitCountsPerOrg();
    my %nTrimPerOrg;
    my $sumCount = 0;
    foreach my $org (sort keys %BaitCountsPerOrg){
        $nTrimPerOrg{$org}=0;
        my $count = $BaitCountsPerOrg{$org}->{max};
        $sumCount += $count;
        $nTrimPerOrg{$org} = $count - $maxBaitsPerOrg if($count > $maxBaitsPerOrg);
    }
    TrimBaitsFromOrganisms('Max',($bQuiet) ? $sumCount: undef,%nTrimPerOrg);
    my $PreMaxBaits = $AbsMaxBaits;
    $AbsMaxBaits = GetBaitCount('Max',(values %RegionOrder));
    $Logger->Log("Filtered out @{[$PreMaxBaits - $AbsMaxBaits]} Baits",$logLevel);
}

sub TrimBaits($$){
    $Logger->Log("TrimBaits","DEBUG");
    my ($densMode,$nBaitsToTrim) = (@_);
    my $TotalTrimmed = 0;
    my $bCycle = 0;
    my $remaining = $nBaitsToTrim;
    while(!$bCycle and $remaining > 0){ #Repeat until the desired number of baits have been removed or no more baits can be removed
        $bCycle = 1;
        my %BaitCountsPerOrg = GetBaitCountsPerOrg();
        $BaitCountsPerOrg{$_} = $BaitCountsPerOrg{$_}->{lc($densMode)} foreach(sort keys %BaitCountsPerOrg);
        my @Orgs = sort {$BaitCountsPerOrg{$b} <=> $BaitCountsPerOrg{$a}} (keys %BaitCountsPerOrg);
        my %nTrimPerOrg;
        @nTrimPerOrg{@Orgs} = (0) x @Orgs;
        for(my $i = 0; $i < @Orgs; $i++){
            my $modulus = 0;
            my $diff = 0;
            if($i < $#Orgs and $BaitCountsPerOrg{$Orgs[$i+1]} > $MIN_BAITS_PER_ORG){
                $diff = $BaitCountsPerOrg{$Orgs[$i]} - $BaitCountsPerOrg{$Orgs[$i + 1]};
            } else {
                $diff = $BaitCountsPerOrg{$Orgs[$i]} - $MIN_BAITS_PER_ORG;
            }
            $diff = 0 if($diff < 0);
            if($diff * ($i + 1) > $remaining){
                $modulus = $remaining % ($i+1);
                $diff = int($remaining/($i+1));
            }
            for(my $j = 0; $j <= $i; $j++){
                my $maxTrim = $BaitCountsPerOrg{$Orgs[$j]} - $MIN_BAITS_PER_ORG;
                $maxTrim = ($maxTrim > 0) ? $maxTrim : 0;
                my $availableTrim = ($maxTrim > $nTrimPerOrg{$Orgs[$j]}) ? $maxTrim - $nTrimPerOrg{$Orgs[$j]} :0;
                my $ToTrim = ($availableTrim > $diff + $modulus) ? ($diff + $modulus) : $availableTrim;
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
        foreach my $org (sort keys %nTrimPerOrg){
            $Val{sum} += $nTrimPerOrg{$org};
            $Val{max} = $nTrimPerOrg{$org} if($Val{max} < $nTrimPerOrg{$org});
            $Val{min} = $nTrimPerOrg{$org} if($Val{min} > $nTrimPerOrg{$org});
        }
        my $trimmed = TrimBaitsFromOrganisms($densMode,$nBaitsToTrim,%nTrimPerOrg);
        $TotalTrimmed += $trimmed;
        $remaining = $nBaitsToTrim - $TotalTrimmed;
        $bCycle = 1 if($trimmed == 0);
    }
    if($TotalTrimmed < $nBaitsToTrim){
        $Logger->Log("Could not trim $nBaitsToTrim baits due to region sharing or size: Trimmed $TotalTrimmed baits","WARNING");
    }
}

sub TrimBaitsFromOrganisms($$%){
    #Entire Regions may be trimmed away
    #As The first 3 sorting priorities are unaffected by region length, and ties past that point are
    #unlikely, Baits are successively removed from the lowest priority regions
    $Logger->Log("TrimBaitsFromOrganism","DEBUG");
    my $densMode = shift(@_);
    my $trimCap = shift(@_);
    my %nBaitsPerOrg = @_;
    my %Record = @_;
    my $bTrimCap = 1;
    unless(defined $trimCap){
        $bTrimCap = 0;
        $trimCap = 0;
        $trimCap += $nBaitsPerOrg{$_} foreach(sort keys %nBaitsPerOrg);
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
            my $message = "Could not trim all $Record{$org} baits from $org due to excessive sharing:".
                " Trimmed @{[$Record{$org} - $nBaitsPerOrg{$org}]}";
                $Logger->Log($message,"WARNING");
        }
    }
    return $TotalTrimmed;
}

sub RemoveRegion($){
    $Logger->Log("RemoveRegion","DEBUG");
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
    $Logger->Log("TrimBaitsFromRegion","DEBUG");
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
    $Logger->Log("Optimizing max baits per organism...","INFO");
    my $nMinClass = 2;
    my @MinOrgs;
    for(my $i = 0; $i < $nMinClass; $i++){
        $nMinClass = 0;
        my (undef,%OrgTilingClass) = DetermineTilingClass();
        my %BaitCountPerOrg = GetBaitCountsPerOrg();
        $BaitCountPerOrg{$_} = $BaitCountPerOrg{$_}->{max} foreach(sort keys %BaitCountPerOrg);
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
        $Logger->Log("Optimization Iteration @{[$i+1]}: Max @ $maxBaitsPerOrg","DEBUG");
    }
    $Logger->Log("Found optimum at $maxBaitsPerOrg baits per organism","INFO");
}

sub OptimizeMaxRegionsPerOrg(){
    $Logger->Log("Optimizing max regions per organism...","INFO");
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
                $maxRegPerOrg = scalar(@{$RegionOrder{$MinOrgs[$i]}});
            }
            ApplyMaxRegionFilter(%RegionOrder,1);
        } elsif(!$maxRegPerOrg){
            $maxRegPerOrg = $maxBaits;
        }
        $Logger->Log("Optimization Iteration @{[$i+1]}: Max @ $maxRegPerOrg","DEBUG");
    }
    $maxRegPerOrg = 1 if($maxRegPerOrg < 1);
    $Logger->Log("Found optimum at $maxRegPerOrg regions per organism","INFO");
}

sub ParseConfig($){
    my $file = shift;
    $file = "$FindBin::RealBin/../HUBDesign.cfg" unless defined $file;
    unless(-e $file){
        warn "Could not find Config file ($file): Revert to Hardcoded defaults\n";
        return;
    }
    my %Config = LoadConfig($file,"WARNING");
    my $_ID_PADDING = $Config{'ID-padding'} if(exists $Config{'ID-padding'}); 
    my $MIN_BAITS_PER_ORG = $Config{'baits-per-org-min'} if(exists $Config{'baits-per-org-min'}); 
    $DEFAULT{d} = $Config{'tiling-density'} if(exists $Config{'tiling-density'}); 
    $DEFAULT{l} = $Config{'length'} if(exists $Config{'length'}); 
}

