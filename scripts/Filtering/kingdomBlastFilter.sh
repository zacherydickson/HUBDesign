#!/bin/bash

if [ "$#" -lt 5 ]; then 
    >&2 echo "Usage: $(basename $0) Threads reshapedOligFile DupInfo BlastDB ValidKingdom [ExcludeList] > ExclusionIntervals.tab";
    exit 1;
fi

Threads=$1;
OligFile=$2;
DupInfo=$3;
BlastDB=$4;
ValidKingdom="$5";
ExcludeList=$6;
tmpEx="/tmp/ex_$(rndstr 8)";
if [ -z "$ExcludeList" ];then
    touch "$tmpEx"
    ExcludeList="$tmpEx"
fi

FastaFile="/tmp/fasta_$(rndstr 8)";
>&2 echo "FastaFile: $FastaFile";

awk -F '\\t' '
    (ARGIND == 1){Exclude[$1]=1;next}
    (ARGIND == 2){Dup[$1]=$2;next}
    (Seen[Dup[$1]] == 1) {next}
    (Dup[$1]){
        Seen[Dup[$1]] = 1; 
        $1 = Dup[$1];
    }
    (!Exclude[$1]){
        print ">"$1;
        print $4;
    }
' "$ExcludeList" "$DupInfo" "$OligFile" > "$FastaFile";
if [ "$?" -ne 0 ]; then 
    >&2 echo "ERROR in creating Fasta File\n";
    rm -f "$tmpEx";
    exit 1;
fi
rm -f "$tmpEx";

blastn -query "$FastaFile" -task blastn -db "$BlastDB" -evalue 0.01 -outfmt '6 qseqid length sskingdom qstart qend' \
    -perc_identity 75 -num_threads "$Threads" | 
    awk -v kingdom="$ValidKingdom" '($2 < 30 || $3 == kingdom || $3 == "N/A") {next} {print $1"\t"$4"\t"$5;}'
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR Blast Searching or Parsing\n";
    exit 1;
fi

rm -f "$FastaFile";
