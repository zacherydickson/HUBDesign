package HUBDesign::Grouping;
use strict;
use Carp;

#Data structure describing groupings seen in gene clusters or a guide tree
#Has 3 public members:
#   members - The ids of all taxa in the grouping
#   count   - The number of times the grouping has been observed
#   isObserved - Distinguishes groupings actually observed in the data from those implied
#   by exclusion
#Has 3 public methods:
#   increment - Increases the count by 1 (or an optional step value); returns the new count
#   isCompatible - Tests if two groupings are compatible; returns boolean
#   isSubset - Tests wether the first grouping is a subset of the latter; returns boolean

sub new {
    my ($package, %args) = @_;
    my $members = (exists $args{members}) ? $args{members} : [];
    my $count = (exists $args{count}) ? $args{count} : 0;
    my $bObs = (exists $args{bObs}) ? $args{bObs} : 1;
    my $self = {'Grouping::Mem' =>$members, 'Grouping::Count' => $count, 'Grouping::Obs' => $bObs};
    bless $self =>$package;
    return $self;
}

sub members {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Grouping'){
        carp "Attempt to call Grouping::members with a non-HUBDesign::Grouping object\n";
    }
    my $memRef = shift;
    if(defined $memRef){
        unless(ref $memRef eq 'ARRAY'){
            croak "Attempt to alter the membership of a HUBDesign::Grouping Object with a non-ARRAY reference\n";
        }
        $self->{'Grouping::Mem'} = $memRef;
    }
    carp "Too many arguments to Grouping::members\n" if(@_);
    return @{$self->{'Grouping::Mem'}};
}

sub count {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Grouping'){
        carp "Attempt to call Grouping::count with a non-HUBDesign::Grouping object\n";
    }
    my $count = shift;
    if(defined $count){
        $self->{'Grouping::Count'} = $count;
    }
    carp "Too many arguments to Grouping::count\n" if(@_);
    return $self->{'Grouping::Count'};
}

sub isObserved {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Grouping'){
        carp "Attempt to call Grouping::isObserved with a non-HUBDesign::Grouping object\n";
    }
    my $bObs = shift;
    if(defined $bObs){
        $self->{'Grouping::Obs'} = $bObs;
    }
    carp "Too many arguments to Grouping::isObserved\n" if(@_);
    return $self->{'Grouping::Obs'};
}

sub increment {
    my $self = shift;
    unless(defined $self and ref $self eq 'HUBDesign::Grouping'){
        carp "Attempt to call Grouping::increment with a non-HUBDesign::Grouping object\n";
    }
    my $step = shift;
    carp "Too many arguments to Grouping::increment\n" if(@_);
    $step = 1 unless(defined $step);
    $self->{'Grouping::Count'} += $step;
    return $self->{'Grouping::Count'}
}

#Tests whether two Grouping objects are under the assumption that a grouping is compatible if
#   The two groupings have no intersection
#   One grouping is a subset of the other
sub isCompatible {
    my $self = shift;
    my $other = shift;
    carp "Too many arguments to Grouping::isCompatible\n" if(@_);
    my $bCompat = 0;
    unless(defined $self and defined $other){
        croak "Attempt to call Grouping::isCompatible with undefined groupings\n";
    }
    unless(ref $self eq 'HUBDesign::Grouping' and ref $other eq 'HUBDesign::Grouping'){
        croak "Attempt to call Grouping::isCompatible with non-HUBDesign::Grouping objects\n";
    }
    my @minorMem = @{$self->{'Grouping::Mem'}};
    my @majorMem = @{$other->{'Grouping::Mem'}};
    if(scalar(@{$self->{'Grouping::Mem'}}) > scalar(@{$other->{'Grouping::Mem'}})){
        my @tmp = @minorMem;
        @minorMem = @majorMem;
        @majorMem = @tmp;
    }
    my %majorSet;
    @majorSet{@majorMem} = (1) x @majorMem;
    my $intersect = 0;
    foreach (@minorMem){
        $intersect++ if(exists $majorSet{$_});
    }
    $bCompat = 1 if(!$intersect or $intersect == scalar(@minorMem));
    return $bCompat;
}

#self is expected to be the smaller of two sets
sub isSubset {
    my $self = shift;
    my $other = shift;
    carp "Too many arguments to Grouping::isSubset\n" if(@_);
    my $bSubset=0;
    unless(defined $self and defined $other){
        croak "Attempt to call Grouping::isSubset with undefined groupings\n";
    }
    unless(ref $self eq 'HUBDesign::Grouping' and ref $other eq 'HUBDesign::Grouping'){
        croak "Attempt to call Grouping::isSubset with non-HUBDesign::Grouping objects\n";
    }
    my @minorMem = @{$self->{'Grouping::Mem'}};
    my @majorMem = @{$other->{'Grouping::Mem'}};
    my %majorSet;
    @majorSet{@majorMem} = (1) x @majorMem;
    my %Union;
    @Union{(@minorMem,@majorMem)} = (1) x (@minorMem + @majorMem);
    $bSubset = 1 if(scalar(keys %Union) == scalar(keys %majorSet));
    return $bSubset;
}

1;
