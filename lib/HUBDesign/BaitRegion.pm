package HUBDesign::BaitRegion;
use strict;
use Carp;


#Data structure representing a bait region
#Has 5 public members: uid, taxon_id, clust_id, pos, and len
#Has the methods: subbr, exclude, toStr, isContig,
#                 isOverlapping, merge, copy

sub new(){
    my $package = shift;
    my $self = {};
    my %args = @_;
    $self->{'BaitRegion::uid'} = $args{uid} if(exists $args{uid});
    delete $args{uid};
    $self->{'BaitRegion::taxon_id'} = $args{taxon_id} if(exists $args{taxon_id});
    delete $args{taxon_id};
    $self->{'BaitRegion::clust_id'} = $args{clust_id} if(exists $args{clust_id});
    delete $args{clust_id};
    $self->{'BaitRegion::pos'} = $args{pos} if(exists $args{pos});
    delete $args{pos};
    $self->{'BaitRegion::len'} = $args{len} if(exists $args{len});
    delete $args{len};
    if(exists $args{seq}){
        $self->{'BaitRegion::seq'} = $args{seq};
        $self->{'BaitRegion::len'} = length($args{seq});
        delete $args{seq};
    }
    carp("Unknown parameters to BaitRegion::new") if(scalar(keys %args));
    bless $self => $package;
    return $self;
}

sub uid(){
    my $self = shift;
    if(@_){
        $self->{'BaitRegion::uid'} = shift;
        carp("Too many arguments to BaitRegion::uid") if (@_);
    }
    return $self->{'BaitRegion::uid'};
}

sub taxon_id(){
    my $self = shift;
    if(@_){
        $self->{'BaitRegion::taxon_id'} = shift;
        carp("Too many arguments to BaitRegion::taxon_id") if (@_);
    }
    return $self->{'BaitRegion::taxon_id'};
}

sub clust_id(){
    my $self = shift;
    if(@_){
        $self->{'BaitRegion::clust_id'} = shift;
        carp("Too many arguments to BaitRegion::clust_id") if (@_);
    }
    return $self->{'BaitRegion::clust_id'};
}

sub pos(){
    my $self = shift;
    if(@_){
        $self->{'BaitRegion::pos'} = shift;
        carp("Too many arguments to BaitRegion::pos") if (@_);
    }
    return $self->{'BaitRegion::pos'};
}

sub seq(){
    my $self = shift;
    if(@_){
        $self->{'BaitRegion::seq'} = shift;
        $self->{'BaitRegion::len'} = length($self->{'BaitRegion::seq'});
        carp("Too many arguments to BaitRegion::seq") if (@_);
    }
    return $self->{'BaitRegion::seq'};
}

sub len(){
    my $self = shift;
    if(@_){
        if(defined $self->{'BaitRegion::seq'}){
            carp("BaitRegion::len can only be set if the region's sequence is undefined");
            return $self->{'BaitRegion::len'};
        }
        $self->{'BaitRegion::len'} = shift;
        carp("Too many arguments to BaitRegion::len") if (@_);
    }
    return $self->{'BaitRegion::len'};
}


#Returns a sub Bait region which is a copy of its parent bait region but whose sequence is a substring
#of its parents
sub subbr(){
    my ($self,$start,$end,$uid) = @_;
    unless(defined $start and defined $end and $start <= $end and $start >= 1 and $end <= length($self->seq)){
        croak "Start and End must be in order and within the sequence in BaitRegion::subbr\n";
    }
    $uid = $self->uid unless(defined $uid);
    return HUBDesign::BaitRegion->new(uid => $uid, taxon_id => $self->taxon_id, clust_id => $self->clust_id, pos => $self->pos + $start - 1, seq => substr($self->seq,$start-1,$end - $start + 1));
}

#Given an reference to list of intervals to exclude from a bait region, returns a list of 
#new bait regions formed by breaking up the original
#Assumes the input intervals are non-overlapping
sub exclude(){
    my ($self,$ivListRef) = @_;
    my @children;
    @{$ivListRef} = sort {$a->[0] <=> $b->[0]} @{$ivListRef};
    my $offset = 0;
    my $origLen = length($self->seq);
    foreach my $iv (@{$ivListRef}){
        my ($start, $end) = map {$_ - $offset} @{$iv};
        $end--;
        unless(defined $start and defined $end and $start <= $end and $start >= 1 and $end <= length($self->seq)){
            croak "Start($start) and End($end) must be in order and within the sequence(1-".length($self->seq).") in BaitRegion::exclude \n";
        }
        if($start > 1){
            push(@children,$self->subbr(1,$start -1,$self->uid."_".join(":",($offset+1,$iv->[0]-1))));
        }
        if($end < length($self->seq)){
            $self = $self->subbr($end +1,length($self->seq));
            $offset += $end;
        } else{
            $self = undef;
            last;
        }
    }
    if(defined $self){
        $self->uid($self->uid."_".join(":",($offset+1,$origLen)));
        push(@children,$self);
    }
    return @children;
}


#Converts the Bait Region to a string for printing
sub toStr(){
    my $self = shift;
    my $bUID = shift || 0;
    my $first = $self->taxon_id;
    $first = $self->uid."\t".$first if($bUID);
    return join("\t",($first,$self->clust_id,$self->pos,$self->len,$self->seq));
}

#Checks if two bait regions are contiguous with eachother, i.e. the first bp of the second region is
#the second bp of the first region
sub isContig(){
    my ($self,$other) = @_;
    unless(ref($self) eq "HUBDesign::BaitRegion" and ref($other) eq "HUBDesign::BaitRegion"){
        carp "[WARNING] Cannot check contiguity of non BaitRegion objects\n";
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

sub isOverlapping(){
    my ($self,$other) = @_;
    unless(ref($self) eq "HUBDesign::BaitRegion" and ref($other) eq "HUBDesign::BaitRegion"){
        carp "[WARNING] Cannot check overlap of non BaitRegion objects\n";
        return 0;
    }

    if($self->taxon_id eq $other->taxon_id){
        if($self->clust_id eq $other->clust_id){
            ($self,$other) = ($other,$self) if($other->pos < $self->pos);
            if($self->pos + $self->len - 1 >= $other->pos){
                return 1;
            }
        }
    }
    return 0;
}

#returns a new object which is the union of the two objects
sub merge(){
    my ($self,$other) = @_;
    unless(ref($self) eq "HUBDesign::BaitRegion" and ref($other) eq "HUBDesign::BaitRegion"){
        carp "[WARNING] Cannot merge non BaitRegion objects\n";
        return undef;
    }
    unless($self->taxon_id eq $other->taxon_id and 
        $self->clust_id eq $other->clust_id){
        carp "[WARNING] Cannot merge BaitRegions originating from different taxa or clusters";
        return undef;
    }

    ($self,$other) = ($other,$self) if($other->pos < $self->pos);
    my $extend = $other->pos + $other->len - $self->pos - $self->len;
    $extend = 0 if($extend < 0);
    my $obj = HUBDesign::BaitRegion->new(uid => $self->uid, clust_id => $self->clust_id,
        taxon_id => $self->taxon_id, len => $self->len + $extend,
        pos => $self->pos);
    if(defined $self->seq and defined $other->seq){
        my $seq = $self->seq;
        $seq .= substr($other->seq,-$extend);
        $obj->seq($seq);
    }
    return $obj;
}

sub copy(){
    my $self = shift;
    my $copy = HUBDesign::BaitRegion->new(
        uid => $self->uid,
        clust_id => $self->clust_id,
        taxon_id => $self->taxon_id,
        pos => $self->pos,
        len => $self->len
    );
    $copy->seq($self->seq) if(defined $self->seq);
    return $copy;
}

1;
