#!/bin/bash

if [ "$#" -lt 5 ];then
    >&2 echo "Usage: $(basename $0) MinDensity MaxBaits POInfo RegionData Tree PriorityClusters [TilingInfo] > Baits.fna";
    exit 1;
fi

MinDensity=$1;
MaxBaits=$2;
POData=$3;
RegionData=$4;
Tree=$5
PriorityClusters=$6;
infoFile=$7;

tileDir="/home/zac/scripts/Hendrik/IntestinalParasites/Redesign";
scriptDir="/home/zac/scripts/Hendrik/COVID";

WorkingDir="/scratch/zac/tile_$(rndstr 8)"
>&2 echo "WorkingDir: $WorkingDir";
mkdir -p "$WorkingDir";

seqFile="$WorkingDir/fasta";
if [ -z "$infoFile" ]; then
    infoFile="$WorkingDir/info";
fi

awk -F '\\t' -v seqFILE="$seqFile" '
    (ARGIND == 1) {
        #Start[$1";"$2] = $4;
        next;
    }
    {
        #gene=$2;
        txid=$2;
        #i=1;
        #while(Start[sprintf("Cluster%04d",i)";"gene] && Start[sprintf("Cluster%04d",i)";"gene] <= $3){
        #    i++;
        #}
        #cluster=sprintf("Cluster%04d",i-1)
        #if(Start[cluster";"gene]){
            printf("%s\t%s\t%s\t%d\n",$1,txid,"Arb",$3);
            print ">"$1 > seqFILE;
            print $4 > seqFILE;
        #} else {
        #    system(">&2 echo Could not find Cluster for "$1);
        #}
    }
' "$POData" "$RegionData" > "$infoFile";

perl "$scriptDir/createTileInfo.pl" "$POData" "$RegionData" "$Tree" > "$infoFile";
if [ "$?" -ne 0 ]; then >&2 echo "ERROR in creating Tile Info"; exit 1; fi

awk '{print ">"$1; print $4;}' "$RegionData" > "$seqFile"

perl "$tileDir/tileBaits.pl" -p "$PriorityClusters" -d  "$MinDensity" -n "$MaxBaits" -av "$infoFile" "$seqFile"

rm -rf "$WorkingDir";
