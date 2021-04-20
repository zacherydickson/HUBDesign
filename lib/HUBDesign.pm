package HUBDesign;

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(OpenFileHandle GetProcessorCount ProcessNumericOption);
our %EXPORT_TAGS = (All => [qw(&OpenFileHandle &GetProcessorCount &ProcessNumericOption)]);

#Given a file path, opens it for reading and fails with the provided exit function if unsuccessful
# The default exit function is die
#Input	File - Path to the file to open
#	Type - A descriptor of the type of file being opened
# (opt)	exit_Func - a reference to a subroutine to run if opening the file fails
#Output	A file handle if successfule, 0 otherwise.
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

1;
