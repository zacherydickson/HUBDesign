#! /usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Std;
use File::Spec;
use FindBin;
use lib File::Spec->catdir($FindBin::RealBin,'..','lib');
use HUBDesign::Util qw(ValidateThreadCount ProcessNumericOption OpenFileHandle LoadConfig);
use HUBDesign::Logger;

my %opts;
getopts('a:l:h',\%opts);

if(exists $opts{h}){
    die "Usage: ".basename($0). "[-a aligner -l lcr-tool] > Settings.cfg\n".
        "\tGenerates a HUBDesign config file given the provided executable names\n".
        "\tAttempts to find the blastn executable as well\n".
        "\tSupported Aligners: [mafft], muscle, clustalo\n".
        "\t\tAligners must be able to output to stdout\n".
        "\tSupported LCR filters: [dustmasker], sdust\n".
        "\tAll other config options must be set manually\n";
}

$opts{a} = (exists $opts{a}) ? lc($opts{a}) : 'mafft';
$opts{l} = (exists $opts{l}) ? lc($opts{l}) : 'dustmasker';

my $Logger = HUBDesign::Logger->new(level => "WARNING");

my $algn_exec = `which $opts{a}`;
unless($algn_exec){
    $Logger->Log("Could not locate $opts{a} executable: manually enter path under Aligner-exec in cfg file","WARNING");
    $algn_exec = $opts{a};
}

my $aligner_cmd = "--thread {threads} -in {input}";
if($opts{a} eq 'mafft'){
    $aligner_cmd = "--auto --localpair --quiet --thread {threads} {input}";
} elsif($opts{a} eq 'muscle'){
    $aligner_cmd = "-quiet -in {input}";
} elsif($opts{a} eq 'clustalo'){
    $aligner_cmd = "--threads={threads} --infile={input}";
} else {
    $Logger->Log("Unsupportted aligner specified: manually enter cmd format under Aligner-cmd in cfg file\n\tFields for the number of threads {threads} and input file {input} may be used","WARNING");
}

my $blast_exec = `which blastn`;
unless($blast_exec){
    $Logger->Log("Could not locate blastn executable: manually enter path under BLAST-exec in cfg file","WARNING");
    $blast_exec = "blastn";
}

my $lcr_exec = `which $opts{l}`;
unless($lcr_exec){
    $Logger->Log("Could not locate $opts{l} executable: manually enter path under lcr-exec in cfg file","WARNING");
    $lcr_exec = $opts{l};
}

my $lcr_cmd = "-w {window} -t {threshold} -i {input}";
if($opts{l} eq 'dustmasker'){
    $lcr_cmd = "-outfmt acclist -window {window} -level {threshold} -in {input}";
}elsif ($opts{l} eq 'sdust'){
    $lcr_cmd = "-w {window} -t {threshold} {input}";
} else {
    $Logger->Log("Unsupported lcr tool specified: manually enter cmd format under lcr-cmd in cfg file","WARNING");
}

my $config = 
"#This config file expects inputs in the format KEY\tVALUE\n".
"#Empty lines, lines which do not have exactly two fields and lines with a leading '#' are ignored\n".
"##========================\n".
"#Dependencies\n".
"\n".
"#Aligner-exec is the path to the executable for aligner to use during clustering\n".
"#Aligner-cmd is a string detailing the options for an aligner; {threads} and {input} are replaced by corresponding values during execution; These options are provided so another aligner may be used which specifies inputs differently\n".
"#Aligner-format is the output format from the aligner, supported formats can be found in the documentations for Bio::AlignIO  \n".
"Aligner-exec	$algn_exec\n".
"Aligner-cmd	$aligner_cmd\n".
"# Ex with -exec	muscle: -in {input} -quiet\n".
"Aligner-format	fasta\n".
"\n".
"#BLAST-exec is the path to the executable to use during blast filtering\n".
"BLAST-exec	$blast_exec\n".
"\n".
"#BOND-exec is the path to the executable for BOND used during identification, BOND is packaged with HUBDesign and this option is not expected to be necessary, {bin is replaced with the location of executable during execution}\n".
"#Another tool could be used, but the cmd and output format are tightly integrated with HUBDesign\n".
"\n".
"BOND-exec	{bin}/../SA_BOND/sa_bond\n".
"\n".
"#lcr-exec is the path to executable for lcr filtering\n".
"#lcr-cmd is a string detailing the call for the lcr identifier; {window} {threshold} {input} are replaced with the relevant values during execution; These options are provided so that another lcr-identifcation tool may be used; Note: the output of the lcr tool is expected to be tab delimited with 3 columns: id, start and end; the acclist outfmt of dustmasker will work ; Note: spaces are required between flags and values in the specification\n".
"lcr-exec	$lcr_exec\n".
"lcr-cmd	$lcr_cmd\n".
"#Ex with -exec	sdust: -w {window} -t {threshold} {input}\n".
"\n".
"###\n".
"#==========================\n".
"##Default Option Values\n".
"length	75\n".
"r2t-divergence	0.15\n".
"translation-table	1\n".
"weight-guide-tree	1\n".
"penetrance	50\n".
"similarity	75\n".
"contiguity	15\n".
"#A value of zero for candidate-pool indicates that all possible candidate probes are identified during the identification phase\n".
"candidate-pool	0\n".
"gc-min	30\n".
"gc-max	70\n".
"#There are two different ways of specifying melting temperature, explicitly with tm-min and tm-max\n".
"#   or finiding the best interval of a given size with tm-range\n".
"#For configuration only one method should be given a default value, and the other value(s) should be None\n".
"tm-min	None\n".
"tm-max	None\n".
"tm-range	10\n".
"lcr-window	64\n".
"lcr-threshold	20\n".
"evalue	10\n".
"hsp-identity	75\n".
"hsp-length	20\n".
"tiling-density	10\n".
"#A value of zero for threads indicates an attempt to use the system maximum threads will be made, failing that only 1 thread will be used\n".
"threads	0\n".
"\n".
"\n".
"#=========================\n".
"#Non-Command Line options\n".
"#Sequences for alignment are temporarily stored here\n".
"cluster-dir	/tmp\n".
"#The width of the numeric portion of any ids generated\n".
"ID-padding	6\n".
"#The minimum number of baits per org to shoot for during selection\n".
"baits-per-org-min	30\n".
"#When no probe-count is provided the number of genomes is multiplied by this value to determine the default upper limit on probe- set size\n".
"genome-mult	100\n";

print $config;
