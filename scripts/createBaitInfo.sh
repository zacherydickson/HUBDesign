#!/bin/bash

if [ "$#" -lt 2 ]; then
    >&2 echo -e "Usage: $(basename $0) TileInfo Sequence\n" \
        "\tTile Info has columns of RegID, txidStr, clustID, clustPos, txid, POPos";
    exit 1;
fi

awk '
    (ARGIND == 1) {
        txid[$1]=$5;
        clustid[$1]=$3;
        POpos[$1]=$6;
        clustpos[$1]=$4;
        next
    }
    /^>/ {
        seqcount++;
        split(substr($0,2),bait,"+");
        printf("Bait%06d\t%s\t%s\t%d\t%s\t%s\n",seqcount,txid[bait[1]],bait[1],POpos[bait[1]]+bait[2],clustid[bait[1]],clustpos[bait[1]]+bait[2])
    }' "$1" "$2"
