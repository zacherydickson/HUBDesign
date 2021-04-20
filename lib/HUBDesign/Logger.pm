package Logger;
use strict;

use Carp;
use Scalar::Util qw(openhandle);

our %_LOG_LEVEL = (ERROR => 0, WARNING => 1, INFO => 2, DEBUG => 3);
our $_DEFAULT_LEVEL = "WARNING";

sub new {
    my ($package, %Args) = @_;
    my $loglevel = (exists $Args{level}) ? uc($Args{level}) : $_DEFAULT_LEVEL;
    my $handle = (exists $Args{handle}) ? $Args{handle} : \*STDERR;
    unless(exists $_LOG_LEVEL{$loglevel}){
        my @keys = sort keys %_LOG_LEVEL;
        croak("Logger level must be in @keys");
    }
    delete $Args{level};
    delete $Args{handle};
    carp("Unrecognized arguments to Logger::new") if(scalar(keys %Args));

    my $self = {'Logger::level' => $_LOG_LEVEL{$loglevel}, 'Logger::handle' => $handle};
    bless $self => $package;

    return $self;
}

#Gets/sets the level of logs which will be written
sub level (){
    my $self = shift;
    my $level;
    if(@_){
        $level = uc(shift);
        carp("Too many arguments to Logger::level") if(@_);
    }
    if(defined $level){
        unless(exists $_LOG_LEVEL{$level}){
            my @keys = sort keys %_LOG_LEVEL;
            croak("Logger level must be in @keys");
        }
        $self->{'Logger::level'} = $_LOG_LEVEL{$level};
    }
    return $self->{'Logger::level'};
}

#Gets/sets the handle to which logs should be written
sub handle(){
    my $self = shift;
    my $handle;
    if(@_){
        $handle = openhandle(shift);
        unless(defined $handle){
            croak("Logger::handle must be provided an open file handle");
        }
        carp("Too many arguments to Logger::handle") if(@_);
    }
    if(defined $handle){
        $self->{'Logger::handle'} = $handle;
    }
    return $self->{'Logger::level'};
}

#Writes a message to STDERR with a timestamp and the level of severity, exits if necessary
sub Log(){
    my ($self, $message, $level) = @_;
    croak("Logger::Log requires a message") unless defined $message;
    $level = (defined $level) ? uc($level) : $_DEFAULT_LEVEL;
    carp("Too many arguments to Logger:Log") if(@_ > 3);
    return if($_LOG_LEVEL{$level} > $self->{'Logger::level'});
    my ($sec,$min,$hour) = localtime;
    my $prevHandle = select($self->{'Logger::handle'});
    printf("%02d:%02d:%02d - [%s] %s\n",$hour,$min,$sec,$level,$message);
    exit 1 if($_LOG_LEVEL{$level} == 0);
    select($prevHandle);
}

#Given a hash of parameters, outputs a summary
sub LogParameters(){
    my ($self,%params) = @_;
    return if($self->{'Logger::level'} < $_LOG_LEVEL{INFO});
    my $prevHandle = select($self->{'Logger::handle'});
    print join("",(('=') x 60)),"\n";
    print "Initialized on ".(localtime)."\n";
    print "Parameters:\n";
    my $nchar = 10;
    foreach (keys %params){
        $nchar = length($_) if(length($_) > $nchar);
    }
    foreach my $param (sort keys %params){
        printf("%*s : %s\n",$nchar + 1,$param,$params{$param});
    }
    print join("",(('=') x 60)),"\n";
    select($prevHandle);
}

1; #Exits properly
