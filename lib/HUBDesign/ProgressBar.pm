#!/usr/bin/env perl
package ProgressBar;
use strict;
use Carp;

our @_PROGRESS_LEVELS = qw(_ \ | / - =);

sub new {
    my ($package, %args) = @_;
    my $state = (exists $args{min}) ? $args{min} : 0;
    my $max = (exists $args{max}) ? $args{max} : 100;
    my $name = (exists $args{title}) ? $args{title} : "Progress";
    my @chars = (' ') x 100;
    my $self = {'PB::min' => $state, 'PB::state' => $state, 'PB::max' => $max, 'PB::chars' => \@chars, 'PB::title' => $name};
    bless $self => $package;
    return $self;
}

sub Update {
    my ($self, $val) = @_;
    unless(defined $self->{'PB::chars'}){
        croak "Cannot Update closed progressbar";
    }
    my $curState = int(($val - $self->{'PB::min'}) / $self->{'PB::max'} * 100 * scalar(@_PROGRESS_LEVELS));
    if($curState > $self->{'PB::state'}){
        my $j = int(($curState - 1) / scalar(@_PROGRESS_LEVELS));
        my $i = int($self->{'PB::state'} / scalar(@_PROGRESS_LEVELS));
        while($i < $j){
            $self->{'PB::chars'}->[$i] = $_PROGRESS_LEVELS[$#_PROGRESS_LEVELS];
            $self->{'PB::state'} += scalar(@_PROGRESS_LEVELS);
            $i++; 
        }
        $self->{'PB::state'} = $curState;
        $self->{'PB::chars'}->[$j] = $_PROGRESS_LEVELS[($curState -1) % scalar(@_PROGRESS_LEVELS)];
        $self->{'PB::chars'}->[$j+1] = ">" if($j < $#{$self->{'PB::chars'}});
        printf STDERR "%s\033[K\n|%s|\033[F",$self->{'PB::title'},join('',@{$self->{'PB::chars'}});
    }
}

sub Close {
    my $self = shift;
    for(my $i = 0; $i < @{$self->{'PB::chars'}}; $i++){
        $self->{'PB::chars'}->[$i] = $_PROGRESS_LEVELS[$#_PROGRESS_LEVELS];
    }
    printf STDERR "%s:\033[K\n|%s|\n",$self->{'PB::title'},join('',@{$self->{'PB::chars'}});
    $self->{'PB::chars'} = undef;
}

return 1;
