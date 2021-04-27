package HUBDesign::BaitRegion;
use strict;
use Class::Struct;
use Carp;


#Data structure rperesenting a bait region
#Has 4 public members
#Has the methods: subbr, exclude, toStr, isContig,
struct (BaitRegion => {taxon_id => '$', clust_id => '$', pos => '$', seq => '$'});


#Returns a sub Bait region which is a copy of its parent bait region but whose sequence is a substring
#of its parents
sub subbr(){
    my ($self,$start,$end) = @_;
    unless(defined $start and defined $end and $start <= $end and $start >= 1 and $end <= length($self->seq)){
        croak "Start and End must be in order and within the sequence in BaitRegion::subbr\n";
    }
    return BaitRegion->new(taxon_id => $self->taxon_id, clust_id => $self->clust_id, pos => $self->pos + $start - 1, seq => substr($self->seq,$start-1,$end - $start + 1));
}

#Given an list of intervals to exlculde form a bait region, returns a list of 
#new bait regions formed by breaking up the original
#Assumes the input intervals are non-overlapping
sub exclude(){
    my ($self,$ivListRef) = @_;
    my @children;
    @{$ivListRef} = sort {$a->[0] <=> $b->[0]} @{$ivListRef};
    my $offset = 0;
    foreach my $iv (@{$ivListRef}){
        my ($start, $end) = map {$_ - $offset} @{$iv};
        $end--;
        unless(defined $start and defined $end and $start <= $end and $start >= 1 and $end <= length($self->seq)){
            croak "Start($start) and End($end) must be in order and within the sequence(1-".length($self->seq).") in BaitRegion::exclude \n";
        }
        if($start > 1){
            push(@children,$self->subbr(1,$start -1));
        }
        if($end < length($self->seq)){
            $self = $self->subbr($end +1,length($self->seq));
            $offset += $end;
        } else{
            $self = undef;
            last;
        }
    }
    push(@children,$self) if(defined $self);
    return @children;
}


#Converts the Bait Region to a string for printing
sub toStr(){
    my $self = shift;
    return join("\t",($self->taxon_id,$self->clust_id,$self->pos,$self->seq));
}

#Checks if two bait regions are contiguous with eachother, i.e. the first bp of the second region is
#the second bp of the first region
sub isContig(){
    my ($self,$other) = @_;
    unless(ref($self) eq "BaitRegion" and ref($other) eq "BaitRegion"){
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

1;
