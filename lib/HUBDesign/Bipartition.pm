package HUBDesign::Bipartition;
use strict;
use Carp;

#Class describing a bipartion of the leaves of a tree
#   Taxa are stored as a bit vector with the taxa in the smaller half of the
#       tree labelled with ones, and the larger half with zeros
#   A bipartitions within a single instance are assumed to be for a tree with
#   the same tips, as such the number of tips must be initialized prior to
#   creating any objects of this class
#Has 3 public members:
#   members - The Indexes of all taxa in smaller half of the bipartition
#   count - The number of times the bipartition was observed
#   size - The number of taxa in the smaller half of the bipartition
#Has 5 public methods:
#   initialize_nTaxa - Sets the number of taxa for all bipartions
#   new - Create a Bipartition object from a reference to a list of indexes
#   get_bits - Returns the string representing the bit vector
#   increment - Increases the count member by 1 or the specified amount
#   hasMember - a boolean test of wether the index specified is in the
#       smaller half
#   isCompatible - a boolean test of wether two bipartitions are compatible
#   isSubset - a boolean test of whether a grouping within the first bipartion is
#       a subset of a grouping within the second; which specific grouping is
#       specified with a 2 bit flag

use constant CHAR_WIDTH => 8;
#Implementing this as a package variable means only bips of one size can be handled at any time; changing this value mid program would invalidate all objects made previously
our $NUMBER_OF_TAXA = undef;

sub initialize_nTaxa {
    my ($package, $nTaxa) = @_;
    if(defined $NUMBER_OF_TAXA){
        carp("The number of taxa for the Bipartions is being reset set");
    }
    $NUMBER_OF_TAXA = $nTaxa;
}

sub new {
    my ($package, %args) = @_;
    unless(defined $NUMBER_OF_TAXA){
        croak("The number of taxa must be initialized with a call to Bipartion::initialize_nTaxa, before calls to Bipartion->new");
    }
    my $count = (exists $args{count}) ? $args{count} : 0;
    my $indexListRef = (exists $args{members}) ? $args{members} : [];
    unless(ref($indexListRef) eq 'ARRAY'){
        croak("Members of a Bipartition object must be an array reference");
    }
    my $bitstring = '';
    my $max = -1;
    my $grpSize = scalar(@{$indexListRef});
    #Build the bitvector
    foreach my $index (@{$indexListRef}){
        vec($bitstring,$index,1) = 0b1;
        $max = $index if($index > $max);
    }
    #Ensure the smaller grouping of a bip is always stored
    #In the case of a tie, then the bip with the rightmost 0 is stored
    if($grpSize >= $NUMBER_OF_TAXA / 2){
        @{$indexListRef} = sort {$a <=> $b} @{$indexListRef};
        if($grpSize > $NUMBER_OF_TAXA / 2  || $indexListRef->[$grpSize -1] == $NUMBER_OF_TAXA -1){
            vec($bitstring,$NUMBER_OF_TAXA-1,1) = 0b0 unless($indexListRef->[$grpSize -1] == $NUMBER_OF_TAXA -1);
            $bitstring = ~$bitstring;
            foreach my $i ($NUMBER_OF_TAXA .. (length($bitstring)*CHAR_WIDTH - 1)){
                vec($bitstring,$i,1) = 0b0;
            }
            $max = $NUMBER_OF_TAXA - 1;
            $grpSize = $NUMBER_OF_TAXA - $grpSize;
        }
    }
    my $self = {'Bipartition::Count' => $count, 'Bipartition::Bitstring' => $bitstring, 'Bipartition::MaxIndex' => $max, 'Bipartition::Size' => $grpSize};
    bless $self => $package;
    return $self;
}

sub members {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::members with a non-Bipartition object\n";
    }
    carp "Too many arguments to Bipartition::members\n" if(@_);
    my @indexList;
    for(my $i = 0; $i < length($self->{'Bipartition::Bitstring'}) * CHAR_WIDTH; $i++){
        push(@indexList,$i) if(vec($self->{'Bipartition::Bitstring'},$i,1));
    }
    return @indexList;
}

sub size {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::members with a non-Bipartition object\n";
    }
    carp "Too many arguments to Bipartition::size\n" if(@_);
    return $self->{'Bipartition::Size'};
}

sub count {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::count with a non-Bipartition object\n";
    }
    my $count = shift;
    if(defined $count){
        $self->{'Bipartition::Count'} = $count;
    }
    carp "Too many arguments to Bipartition::count\n" if(@_);
    return $self->{'Bipartition::Count'};
}

sub get_bits {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::get_bits with a non-Bipartition object\n";
    }
    return $self->{'Bipartition::Bitstring'};
}



sub increment {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::increment with a non-Bipartition object\n";
    }
    my $step = shift;
    carp "Too many arguments to Bipartition::increment\n" if(@_);
    $step = 1 unless(defined $step);
    $self->{'Bipartition::Count'} += $step;
    return $self->{'Bipartition::Count'}
}

sub hasMember(){
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::hasMember with a non-Bipartition object\n";
    }
    my $index = shift;
    unless(defined $index){
        croak "Attempt to call Bipartition::hasMember without an index\n";
    }
    my $bNeg = shift;
    $bNeg = 0 unless(defined $bNeg);
    my $bitstring = $self->{'Bipartition::Bitstring'};
    $bitstring = ~$bitstring if($bNeg);
    return vec($bitstring,$index,1);
}

#Two bipartitions A and B are compatible if any of the following are true:
#   There is no intersection between A and B
#   There is no intersection between A and B^c
#   There is no intersection between A^c and B
#   There is no intersection between A^c and B^c (Only required if The on half is not gauranteed to be
#       the smaller half of the bip)
sub isCompatible {
    my $self = shift;
    my $other = shift;
    carp "Too many arguments to Bipartition::isCompatible\n" if(@_);
    my $bSubset=0;
    unless(defined $self and defined $other){
        croak "Attempt to call Bipartition::isCompatible with undefined bipartions\n";
    }
    unless(ref $self eq 'HUBDesign::Bipartition' and ref $other eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::isCompatible with non-Bipartition objects\n";
    }
    #bips with only one (or no) taxa on one side are always compatable with every other bip
    if($self->{'Bipartition::Size'} <= 1 || $other->{'Bipartition::Size'} <= 1){
        return 1;
    }
    my $bCompat = 0;
    my ($nBits) = sort {$b <=> $a} ($self->{'Bipartition::MaxIndex'}, $other->{'Bipartition::MaxIndex'});
    #Below is a sho-off-y way of saying
    #$bCompat = _has_intersect(self,other)
    #$bCompat = _has_intersect(self,~other) if(!$bCompat)
    #$bCompat = _has_intersect(~self,other) if(!$bCompat)
    for(my $test = 0; $test < 3 && !$bCompat; $test++){
        my @bitstrList = ($self->{'Bipartition::Bitstring'},$other->{'Bipartition::Bitstring'});
        #Negate the appropriate bitstring: test=0b00 neither, test=0b01 other, test=0b10 self
        foreach my $i (0 .. 1){
            $bitstrList[$i] = ~$bitstrList[$i] if($test & 1 << $i);
        }
        $bCompat = !_has_intersect(@bitstrList,$nBits);
    }
    return $bCompat;
}

sub _has_intersect {
    my ($bitStr1,$bitStr2,$nbits) = @_;
    my $bIntersect = 0;
    my $intersect = unpack("B*",$bitStr1 & $bitStr2);
    substr($intersect,-1*CHAR_WIDTH,CHAR_WIDTH-($nbits%CHAR_WIDTH)-1,"");
    $bIntersect = ($intersect =~ m/1/) ? 1 : 0;
    return $bIntersect;
}

sub _toStr {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::toStr with a non-Bipartition object\n";
    }
    return join("\t",($self->size,$self->count,unpack("B*",$self->get_bits)));
}

#Tests whether the flag defined self grouping is a subset of the flag defined grouping of
#the other
#Bips are guaranteed to be stored with the small grouping on, and the big grouping off
#Flag bits are 1s bit is self, 2s bit is other
#flag : 0b11 compare big to big, 0b10 compare big to small
#     : 0b01 compare small to big, 0b00 compare small to small
sub isSubset {
    my $self = shift;
    my $other = shift;
    my $flag = shift;
    carp "Too many arguments to Bipartition::isSubset\n" if(@_);
    my $bSubset=0;
    unless(defined $self and defined $other){
        croak "Attempt to call Bipartition::isSubset with undefined bipartions\n";
    }
    unless(ref $self eq 'HUBDesign::Bipartition' and ref $other eq 'HUBDesign::Bipartition'){
        croak "Attempt to call Bipartition::isSubset with non-Bipartition objects\n";
    }
    unless(defined $flag){
        croak "Attempt to call Bipartition::isSubset with an undefined flag\n";
    }
    #A is a subset of B if the intersect of A and B is A :: (A & B == A)
    my @bitstrList = ($self->{'Bipartition::Bitstring'},$other->{'Bipartition::Bitstring'});
    #Need to make sure that Both bitstrings are the same size for negation to work properly
    if(length($bitstrList[0]) < length($bitstrList[1])){
        vec($bitstrList[0],length($bitstrList[1])*CHAR_WIDTH-1,1) = 0b0;
    } elsif(length($bitstrList[0]) > length($bitstrList[1])){
        vec($bitstrList[1],length($bitstrList[0])*CHAR_WIDTH-1,1) = 0b0;
    }
    #obtain the appropriate grouping through negation
    foreach my $i (0 .. 1){
        if(($flag >> $i) & 1){
            $bitstrList[$i] = ~$bitstrList[$i];
            for(my $j = length($bitstrList[$i])*CHAR_WIDTH; $j > $NUMBER_OF_TAXA; $j--){
                vec ($bitstrList[$i],$j-1,1) = 0b0;
            }
        }
    }
    my $intersect = $bitstrList[0] & $bitstrList[1];
    $bSubset = 1 if($intersect eq $bitstrList[0]); 
    return $bSubset; 
}



return 1;
