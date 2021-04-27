package HUBDesign::Pseudorg;
use strict;
use Scalar::Util qw(looks_like_number);
use Carp;

#Has readonly members: len, seq, count
#Has private members: clustList, posList
#   clustList is a list of CLUSTER Objects
#   posList is the position within the buffered pseudo-organism sequence at which each cluster object
#   begins
#Has functions add, map_pos, new

#Initializes a Pseudorg obj
sub new (){
    my $package = shift;

    my $self = {};
    if(@_){
        my %Args = @_;
        $self->{'Pseudorg::id'} = $Args{id};
        delete($Args{id});
        carp("[WARNING] too many arguments to Pseudorg::init") if(scalar(keys %Args));
    }
    $self->{'Pseudorg::len'} = 0;
    $self->{'Pseudorg::clustList'} = [];
    $self->{'Pseudorg::posList'} = [];
    bless $self => $package;
    return $self;
}

#Gets or sets the id for the object
sub id (){
    my  $self = shift;
    my $id;
    if(@_){
        $id = shift;
        carp("Too many arguments to Pseudorg::id");
    }
    if(defined $id){
        $self->{'Pseudorg::id'} = $id;
    }
    return $self->{'Pseudorg::id'};
}

#Adds a cluster object to the Pseudo-organism
#Returns 0 if unsuccessful, and the current number of clusters otherwise
sub add (){
    my ($self,$clustObj) = @_;
    unless(defined $clustObj){
        carp "[WARNING] Not enough arguments to Pseudorg::add\n";
        return 0;
    }
    unless(ref($clustObj) == "CLUSTER"){
        carp "[WARNING] Cannot Pseudorg::add non-cluster to Pseudorg $self->{'Pseudorg::id'}\n";
        return 0;
    }
    unless(defined $clustObj->uid and defined $clustObj->seq){
        carp "[WARNING] Attempt to add at least partially undefined cluster to Pseudorg $self->{'Pseudorg::id'}\n";
        return 0;
    }
    push(@{$self->{'Pseudorg::clustList'}},$clustObj);
    push(@{$self->{'Pseudorg::posList'}}, $self->{'Pseudorg::len'} + 1);
    $self->{'Pseudorg::len'} += length($clustObj->seq) + 1;
    return scalar(@{$self->{'Pseudorg::clustList'}});
}

#given a 1 indexed position in the buffered pseudo-organism sequence
#returns a two element vector, the id of the cluster that position is in, and the position within that
#cluster
#returns zero if unsuccessful
sub map_pos(){
    my ($self,$pos) = @_;
    unless(defined $pos){
        carp "[WARNING] Not enough arguments to PSEUFORG::map_pos\n";
        return 0;
    }
    unless(looks_like_number($pos) and $pos >= 1 and $pos <= $self->len){
        carp "[WARNING] Cannot Pseudorg::map_pos which is outside the sequence\n";
        return 0;
    }
    my $index = 0;
    $index++ while($index + 1 < $self->count and $self->{'Pseudorg::posList'}->[$index + 1] <= $pos);
    $pos = $pos - $self->{'Pseudorg::posList'}->[$index] + 1;
    #Check if the mapped position happens to the the buffer position between sequences;
    if($pos > length($self->{'Pseudorg::clustList'}->[$index]->seq)){
        $pos = undef;
    } 
    return ($self->{'Pseudorg::clustList'}->[$index]->uid, $pos);
}

#Returns the length of the buffered pseudo-organism sequence
sub len(){
    my $self = shift;
    return $self->{'Pseudorg::len'};
}

#Returns the buffered pseudo-organism sequence
sub seq(){
    my $self = shift;
    return join("n", (map { uc($_->seq) } @{$self->{'Pseudorg::clustList'}}));
}

#Return sthe number of clusters in the psuedo-organism
sub count(){
    my $self = shift;
    return scalar(@{$self->{'Pseudorg::clustList'}});
}

sub get_ids(){
    my $self = shift;
    return map {$_->uid} @{$self->{'Pseudorg::clustList'}};
}

1;
