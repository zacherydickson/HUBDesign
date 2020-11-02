#!/bin/bash

if [ "$#" -lt 2 ]; then
    >&2 echo "Usage: $(basename $0) BONDOutput BONDManifest";
    exit 1;
fi

OligFile=$1
Manifest=$2;

awk -F '\\t' '
    (ARGIND == 1) {ID[FNR]=$0;next}
    {
        baitNo++;
        split($1,ts,":");
        idx=ts[2];
        id = ID[idx];
        split($5,fp,":");
        pos = fp[2]+1;
        print sprintf("Bait%06d",baitNo)"\t"id"\t"pos"\t"$2
    }' "$Manifest" "$OligFile"
