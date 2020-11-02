#!/usr/bin/perl

package DistMat;
use Class::Struct;


struct ( 'DistMat', {lookup => '$', size => '$' , data => '@', label => '@', output => '$'});

sub DistMat::lookup {
    my $self = shift(@_);
    if($_[0] eq 'NEW'){
        return undef;
    }
    my ($i,$j) = @_;
    unless(defined $i and defined $j){
        die "DistMat->lookup requires a row and column\n";
    }
    if($i < 0 or $j < 0){
        die "DistMat->lookup, i and j must be at least 0\n";
    }
    if($i >= $self->size or $j >= $self->size){
        die "DistMat->lookup, i and j must be less than the number of rows\n";
    }
    return 0 if($i == $j);
    ($i,$j) = ($j, $i) if($i > $j);
    my $index = $i*(2*$self->size - $i -3)/2 + $j -1;
    $self->{'DistObj::lookup'} = $self->data->[$index];
    return $self->{'DistObj::lookup'};
}

sub DistMat::output {
    my $self = shift;
    my $fmt = shift;
    $fmt = 'Phylip' unless(defined $fmt);
    if($fmt eq 'NEW'){
        return 0;
    }
    if($fmt eq 'Phylip'){
        print $self->size,"\n";
        for(my $i = 0; $i < $self->size; $i++){
            printf("%-12s",$self->label->[$i]);
            print join(' ',map {$self->lookup($i,$_)} (0 .. ($self->size-1))),"\n";
        }
    }
}

package Main;

use warnings;
use strict;
use File::Basename;
sub LoadMatrix($);

if(@ARGV < 1){
    die "Usage: @{[basename($0)]} distmat\n";
}

my $DistObj = LoadMatrix(shift(@ARGV));
$DistObj->output;

sub LoadMatrix($){
    my $file = shift(@_);
    my @UT;
    my @headers;
    my $bSkip = 1;
    my $bHeader = 1;
    if(open(my $fh, $file)){
        while(my $line = <$fh>){
            chomp($line);
            if($bSkip){
                if($line =~ m/Gap weighting is/){
                    <$fh>;
                    $bSkip = 0;
                }
                next;
            }
            my @values = split(/\s+/,$line);
            shift(@values);
            if($bHeader){
                @headers = @values;
                $bHeader = 0;
                next;
            }
            shift(@values); #Remove Zero from  self comparison
            pop(@values); # Remove ID from line
            pop(@values); # Remove ACC from line
            push(@UT,@values);
        }
        close($fh);
    } else {
        die "Could not open distance matrix: $!\n";
    }
    return DistMat->new(size => scalar(@headers), data => \@UT, label => \@headers, lookup => 'NEW', output => 'NEW');
}



