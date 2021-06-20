#!/usr/bin/perl
package Bio::Align::ParallelDNAStatistics;
use vars qw(%DNAChanges @Nucleotides %NucleotideIndexes
	    $GapChars $SeqCount $DefaultGapPenalty %DistanceMethods
            $CODONS %synchanges $synsites $Precision $GCChhars);
use strict;
use Bio::Align::DNAStatistics;
use Parallel::ForkManager;
use POSIX qw(ceil);
our @ISA = qw(Bio::Align::DNAStatistics);


sub new {
    my ($class,%args) = @_;
    my $self = $class->SUPER::new(%args);
    $self->pairwise_stats( Bio::Align::PairwiseStatistics->new());
    my $threads = $args{'-threads'};
    $threads ||= 1;
    $self->{'Obj::Threads'} = $threads;
    $self->{'Obj::ResCollect'} = {};
    my $pm = Parallel::ForkManager->new(max_proc => $self->{'Obj::Threads'});
    #Each child will return a simple array reference, the identity of each child encodes the covered block
    #of the matrix in the format: i_init,nRow,j_init,nCol
    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $valueListRef) = @_;
            $self->{'Obj::ResCollect'}->{$ident} = $valueListRef;
       }
    );
    $self->{'Obj::ForkManager'} = $pm;
    return $self;
}

#Gets or sets the number of threads for this ParallelDNAStats Object
#It is assumed that the client has validated the number of threads
sub threads {
    my $self = shift;
    if(@_){
        $self->{'Obj::Threads'} = shift;
        if(@_){
            $self->warn("Too many arguments to Bio::Align::ParallelDNAStatistics::threads\n");
        }
    } else {
        return $self->{'Obj::Threads'}
    }
}

sub D_Uncorrected {
    my ($self,$aln,$gappenalty) = @_;
    $gappenalty = $Bio::Align::DNAStatistics::DefaultGapPenalty unless defined $gappenalty;
    return 0 unless $self->_check_arg($aln);
    # ambiguities ignored at this point
    my (@seqs,@names,@values,%dist);
    my $seqct = 0;
    foreach my $seq ( $aln->each_seq) {
        push @names, $seq->display_id;
        push @seqs, uc $seq->seq();
        $seqct++;
    }
    my $precisionstr = "%.$Bio::Align::DNAStatistics::Precision"."f";

    my $len = $aln->length;
     
    #Handle Diagonals where dist is always 0
    for(my $i = 0; $i < $seqct; $i++){
        $dist{$names[$i]}->{$names[$i]} = [$i,$i];
        $values[$i][$i] = sprintf($precisionstr,0);
    }

    #Process Upper triangle calculations in blocks of approximately the same size
    #Blocks are created by selecting the highest, leftmost uncovered element in the UC,
    #then selecting a block which either hits both the right and bottom bounds, or covers
    #at least the desired area.
    my $pm = $self->{'Obj::ForkManager'};
    #Set the number of blocks to be the minimum
    my $nBlocks = (sort {$a <=> $a} ($seqct*($seqct-1)/2,$self->{'Obj::Threads'}))[0];
    my ($top,$left,$right) = (0, 1, $seqct-1);
    my $target = 0.5*$seqct*($seqct-1)/$nBlocks;
    my @coverage;
    for(my $i = 0; $i < $seqct - 1; $i++){
        push(@coverage,[(0)x($seqct)]);
    }
    BLOCK:
    while($top < $seqct -1){
       my $width = $right - $left + 1; 
       my $height = $right - $top;
       #Set maximal block width
       my $nCol = ($width >= $target) ? ceil($target) : $width;
       my $nRow = 1;
       my $adjust = 0;
       #Set minimal block height
       until($nRow >= $height or $nCol * $nRow -$adjust >= $target){
           $nRow++;
           $adjust = $nRow - ($left - $top);
           $adjust = ($adjust < 1) ? 0 : $adjust * ($adjust + 1) / 2;
       }
       #Adjust to minimal block width
       while($nCol * $nRow - $adjust - $nRow >= $target){
           $nCol--;
       }
       for(my $i = $top; $i < $top + $nRow; $i++){
           for(my $j = $left; $j < $left + $nCol; $j++){
               $coverage[$i]->[$j] = 1;
           }
       }
       my $ident = "$top,$nRow,$left,$nCol";
       my $i_init = $top;
       my $j_init = $left;
       while($top < $seqct - 1 and $coverage[$top]->[$left]){
           $left++;
           if($left eq $seqct){
               $top++;
               $left = $top+1;
           }
       }
       $right = $left;
       $right++ until($right >= $seqct -1 or $coverage[$top]->[$right + 1]);
       $pm->start($ident) and next BLOCK;
       my @valueList;
       for(my $i = $i_init; $i < $i_init + $nRow; $i++){
           for(my $j = $j_init; $j < $j_init + $nCol; $j++){
               next if ($i >= $j);
               my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i], $seqs[$j]);
               my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + $matrix->[2]->[2] + $matrix->[3]->[3] );
               my $denom = ( $len - $gaps + ( $gaps * $gappenalty));
               $self->warn("No distance calculated between $names[$i] and $names[$j], default -1") unless $denom;
               my $D = $denom ? 1 - ( $m / $denom) : -1;
               push(@valueList,$D);
           }
       }
       $pm->finish(0,\@valueList);
    }
    $pm->wait_all_children;

    while(my ($ident,$valueListRef) = each %{$self->{'Obj::ResCollect'}}){
        my ($i_init,$nRow,$j_init,$nCol) = split(/,/,$ident);
        my $i_term = $i_init+$nRow -1;
        my $j_term = $j_init+$nCol-1;
        if(defined $valueListRef){
            for(my $i = $i_init; $i <= $i_term; $i++){
                for(my $j = $j_init; $j <= $j_term; $j++){
                    next if($i >= $j);
                    my $D = shift(@{$valueListRef});
                    $dist{$names[$i]}->{$names[$j]} = [$i,$j];
                    $dist{$names[$j]}->{$names[$i]} = [$i,$j];
                    $D = ($D != -1) ? sprintf($precisionstr,$D) : sprintf("%-*s", $Bio::Align::DNAStatistics::Precision + 2, $D);
                    $values[$j][$i] = $values[$i][$j] = $D;
                }
           }
        
        } else {
           $self->warn("No distances calculated for sequences in ($i_init,$j_init) to ($i_term,$j_term)"); 
        }
    }

    $self->{'Obj::ResCollect'} = {};
    return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 
}
