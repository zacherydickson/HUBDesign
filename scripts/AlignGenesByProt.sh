#!/bin/bash

if [ "$#" -lt 3 ]; then
    >&2 echo "Usage: $(basename $0) Threads GenCode DNAFile";
fi

Threads=$1;
GenCode=$2;
DNAFile=$3;

scriptDir="/home/zac/scripts/Hendrik/COVID";

export TMPDIR="/scratch/zac";
#transeq -auto -sequence "$DNAFile" -outseq /dev/stdout -table "$GenCode" -trim 2> /dev/null |
perl ~/scripts/fasta-translate.pl -c "$GenCode" "$DNAFile" |
    mafft --thread "$Threads" --localpair --quiet /dev/stdin |
    perl "$scriptDir/SubCodonForAA.pl" "$DNAFile" /dev/stdin;
