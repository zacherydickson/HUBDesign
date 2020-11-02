#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::Tools::GFF;
use Class::Struct;
use Scalar::Util qw(looks_like_number);

struct ('GeneDat' => {chr => '$', start => '$', end => '$', strand =>'$', ID => '$', version => '$'});

sub GeneDat2str($;$);

if (@ARGV < 1){
    die "Usage: ".basename($0)." gffFile1 ...\n";
}

#Hash with Keys of gene symbols, and values of hash ref with keys of txids and values of array ref with
#elements of GeneDat
my %GeneInfo;

while(my $file = shift(@ARGV)){
    my $txid = basename($file,('.fna.gff','.fasta.gff','.gff'));
    my $IN = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);
    while (my $feature = $IN->next_feature){
        next unless($feature->primary_tag eq 'CDS');
        my $seqID = $feature->seq_id;
        my $start = $feature->start;
        my $end = $feature->end;
        my $strand = $feature->strand;
        my ($ID) = $feature->get_tag_values('ID') if($feature->has_tag('ID'));
        my $gene;
        if ($feature->has_tag('gene')){
            ($gene) = $feature->get_tag_values('gene') 
        } elsif($feature->has_tag('product')){
            ($gene) = $feature->get_tag_values('product');
            if($gene =~ /hypothetical/){
                $gene = "HYPOPROT"
            } else {
                my @words = split(/ |-|,/,$gene);
                $gene = "";
                foreach my $word (@words){
                    next if($word eq 'protein');
                    $word = substr($word,0,2) unless(looks_like_number($word));
                    $gene .= $word;
                }
            }
        }
        my $ver;
        next unless(defined $gene);
        ($gene,$ver) = split(/_/,$gene);
        $ver = 1 unless(defined $ver);
        $gene =~ s/\//_/g;
        $GeneInfo{$gene} = {} unless(exists $GeneInfo{$gene});
        $GeneInfo{$gene}->{$txid} = [] unless(exists $GeneInfo{$gene}->{$txid});
        push(@{$GeneInfo{$gene}->{$txid}},GeneDat->new(chr => $seqID, start => $start,
            end => $end, strand => $strand, ID => $ID, version => $ver));
    }
    $IN->close();
}

foreach my $gene (sort keys %GeneInfo){
    foreach my $txid (sort (keys %{$GeneInfo{$gene}})){
        foreach my $geneObj (@{$GeneInfo{$gene}->{$txid}}){
            print "$gene\t$txid\t".GeneDat2str($geneObj)."\n";
        }
    }
}


sub GeneDat2str($;$){
    my $self = shift;
    my $delim = shift;
    $delim = "\t" unless(defined $delim);
    my @Values = ($self->ID);
    push(@Values,$self->version);
    push(@Values,$self->chr);
    push(@Values,$self->start);
    push(@Values,$self->end);
    push(@Values,$self->strand);
    return(join($delim,@Values));
}
