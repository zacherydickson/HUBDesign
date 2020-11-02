#!/bin/bash

if [ "$#" -lt 1 ]; then
    >&2 echo "Usage: $(basename $0) Root2TipDist ClusterOutput MultipleSequenceAlignment[.gz] > ClusterConsensuses.fna\n";
    exit 1;
fi


scriptDir="/home/zac/scripts/Hendrik/Butyrate";
CollapseMSADir="/home/zac/scripts/Hendrik/Sepsis/BaitDesign";
ReRootDir="/home/zac/scripts/Hendrik/COVID";
WorkingDir="/scratch/zac/Clstr_$(rndstr 8)";

R2TDist=$1;
ClusterFile=$2;
MSAFile=$3;
DistMat="$WorkingDir/distmat";
DistMat="$WorkingDir/distmat";
PhylipMat="$WorkingDir/phylip";
Tree="$WorkingDir/outtree";

genename=${MSAFile%%.*};
genename=${genename##*/};

catCmd="cat";
file "$MSAFile" | grep 'gzip' && catCmd="zcat"

mkdir -p "$WorkingDir";
>&2 echo "WorkingDir: $WorkingDir";

distmat -outfile "$DistMat" -nucmethod 0 <("$catCmd" "$MSAFile") 2>"$WorkingDir/distmat.log";
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in distmat";
    exit 1;
fi

perl "$scriptDir/distmat2phylip.pl" "$DistMat" > "$PhylipMat" 2>"$WorkingDir/dist2phylip.log";
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in conversion to phylip format";
    exit 1;
fi

if [ $(head -1 "$PhylipMat") -gt 2 ]; then
    curDir=$(pwd);
    cd "$WorkingDir"
    neighbor < <(echo -e "phylip\nY\n") > neighbour.log;
    if [ "$?" -ne 0 ]; then
        >&2 echo "ERROR in neighbour create tree";
        exit 1;
    fi
    cd $curDir;
else
    awk '(FNR == 2) {dist=$3; print "(1:"dist/2",2:"dist/2");"}' "$PhylipMat" > "$Tree";
fi

tmpFile=$(rndstr 12);
perl "$ReRootDir/ReRootTree.pl" "$Tree" "$tmpFile"
if [ "$?" -ne 0 ]; then  >&2 echo "Error in Rerooting tree"; exit 1; fi
mv -f "$tmpFile" "$Tree";

perl "$scriptDir/GetClusters.pl" -d "$R2TDist" "$Tree" > "$ClusterFile" 2>"$WorkingDir/GetClusters.log";
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in Get Clusters";
    exit 1;
fi

awk -v wd="$WorkingDir" '
    (ARGIND == 1) {Cluster[$1] = $2;next}
    /^>/ {seqcount++; FILE = wd"/"Cluster[seqcount]".msa"}
    {print $0 >> FILE}' "$ClusterFile" <("$catCmd" "$MSAFile") 2>"$WorkingDir/SplitMSA.log";
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in splitting MSA";
    exit 1;
fi

if [ $(ls "$WorkingDir" | grep -c 'msa$') -lt 1 ]; then
    >&2 echo "ERROR No split MSA Files created";
    exit 1;
fi;

for f in "$WorkingDir"/*.msa; do 
    Cluster=$(basename "$f" ".msa");
    perl "$CollapseMSADir/collapseMSA.pl" -c -a "${genename}_"$(printf "%04d" $Cluster) "$f" 2>"$WorkingDir/CollapseMSA_$Cluster.log";
    if [ "$?" -ne 0 ]; then
        >&2 echo "ERROR in collapse MSA for Cluster $Cluster";
        exit 1;
    fi
done

rm -rf "$WorkingDir";

