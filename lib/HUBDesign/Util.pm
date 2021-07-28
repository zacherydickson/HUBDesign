package HUBDesign::Util;

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(OpenFileHandle GetProcessorCount ProcessNumericOption ValidateThreadCount LoadConfig);
our %EXPORT_TAGS = (All => [qw(&OpenFileHandle &GetProcessorCount &ProcessNumericOption &ValidateThreadCount &LoadConfig)]);

#Given a file path, opens it for reading and fails with the provided exit function if unsuccessful
# The default exit function is die
#Input	File - Path to the file to open
#	Type - A descriptor of the type of file being opened
# (opt)	exit_Func - a reference to a subroutine to run if opening the file fails
#Output	A file handle if successful, 0 otherwise.
#Note: Calling process must close the file handle
sub OpenFileHandle($$;$){
    my ($file,$type,$level) = @_;
    $level = "ERROR" unless(defined $level);
    if(open(my $fh, $file) ){
        return $fh;
    } else {
        my $message = "Could not open $type file ($file): $!";
        my ($sec,$min,$hour) = localtime;
        printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$min,$sec,$level,$message);
        exit 1 if($level eq "ERROR");
        return 0;
    }
}

#Returns the maximum number of processors on a linux based system
#Output: -1 if /proc/cpuinfo not found, the number of processors otherwise
sub GetProcessorCount(){
    my $cpu_count = -1;
    if(open my $handle, "/proc/cpuinfo"){
        $cpu_count = scalar(map /^processor/, <$handle>);
        close($handle);
    }
    return $cpu_count;
}

sub ValidateThreadCount($){
    my $max_proc = shift;
    my $cpu_count = GetProcessorCount();
    my $message = undef;
    if($cpu_count == -1){
        $message = "Could not determine Sys_max threads: $!\n\t Can't fully validate max_threads\n";
        if($max_proc == 0){
            $message = "Sys_max threads unknown: proceeding with 1 thread";
            $max_proc = 1;
        }
    } elsif($max_proc == 0 or $max_proc > $cpu_count){
        if($max_proc > $cpu_count){
            $message = "Max threads is greater than sys_max: proceeding with $cpu_count threads";
        }
        $max_proc = $cpu_count;
    }
    if(defined $message){
        my ($sec,$min,$hour) = localtime;
        printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$min,$sec,"WARNING",$message);
    }
    return $max_proc;
}

sub ProcessNumericOption($$$$$$){
    my($val,$default,$min,$max,$bInt,$varName) = @_;
    return $default unless(defined $val);
    if(looks_like_number($val)){
        $val = int($val) if($bInt);
        if(!defined $min or $val >= $min){
            if(!defined $max or $val <= $max){
                return $val;
            }
        }
    }
    my $message = sprintf("%s must be a%s between %s and %s",$varName,
       ($bInt ? "n integer" : " value"),(defined $min ? $min : "-∞"),(defined $max ? $max : "∞"));
    my ($sec,$minute,$hour) = localtime;
    printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$minute,$sec,"ERROR",$message);
    exit 1;
}

sub LoadConfig($;$){
    my $file = shift;
    my $level = shift;
    $level = "ERROR" unless(defined $level);
    my $fh = OpenFileHandle($file,"Config",$level);
    my %Config;
    while(my $line = <$fh>){
        next if($line =~ /^#/);
        chomp($line);
        next if($line eq "");
        my($key,$value, @other) = split(/\t/,$line);
        next unless(defined $key and defined $value and scalar(@other) == 0);
        if(exists $Config{$key}){
            my $message = "Duplicate config entry for $key on line $. ignored\n";
            my ($sec,$minute,$hour) = localtime;
            printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$minute,$sec,"WARNING",$message);
        }
        $Config{$key} = $value;
    }
    close($fh);
    return %Config;
}

1;
