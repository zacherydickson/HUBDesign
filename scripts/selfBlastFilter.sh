#!/bin/bash

if [ "$#" -lt 3 ]; then 
    >&2 echo "Usage: $(basename $0) Threads reshapedOligFile DupInfo [ExcludeList]";
    exit 1;
fi

Threads=$1;
OligFile=$2;
DupInfo=$3;
ExcludeList=$4;
tmpEx="/scratch/zac/ex_$(rndstr 8)";
if [ -z "$ExcludeList" ];then
    touch "$tmpEx"
    ExcludeList="$tmpEx"
fi

WorkingDir="/scratch/zac/$(rndstr 12)"
mkdir -p "$WorkingDir"
>&2 echo "WorkingDir: $WorkingDir";

FastaFile="$WorkingDir/fasta";
MapFile="$WorkingDir/map";

awk -F '\\t' -v map="$MapFile" '
    BEGIN {lastmap = 0;}
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
        if(!tgt[$2]){lastmap++;tgt[$2]=lastmap;}
        print $1,tgt[$2] >> map;
    }
' "$ExcludeList" "$DupInfo" "$OligFile" > "$FastaFile";
if [ "$?" -ne 0 ]; then 
    >&2 echo "ERROR in creating Fasta File\n";
    rm -f "$tmpEx";
    exit 1;
fi
rm -f "$tmpEx";

BlastDB="$WorkingDir/blast";
makeblastdb -dbtype "nucl" -in "$FastaFile" -out "$BlastDB" -parse_seqids -taxid_map "$MapFile" \
    -blastdb_version 5 > /dev/null;
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in making BlastDB\n";
    exit 1;
fi

Targets=$(tail -1 "$MapFile" | cut -f2 -d' ');
for i in $(seq 1 "$Targets"); do
    awk -v ex="$i" '
        (ARGIND == 1) {Map[$1]=$2;next}
        /^>/ {p=0; if(Map[substr($1,2)] == ex) {p=1}}
        (p==1)
    ' "$MapFile" "$FastaFile" | 
        blastn -task blastn -db "$BlastDB" -evalue 0.01 -outfmt '6 ' \
            -perc_identity 75 -negative_taxids "$i" -num_threads "$Threads";
done #| awk '($2 < 30 || Seen[$1]) {next} {print $1; Seen[$1] = 1;}'

if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR Blast Searching or Parsing\n";
    exit 1;
fi

rm -rf "$WorkingDir";
