#!/bin/bash

if [ "$#" -lt 7 ]; then
    >&2 echo "Usage: $(basename $0) Threads MinPen ClustInfo Tree PsudeoOrgInfoOut ClustAssign ClustSeq1 ... >  PseudoOrgsSeq.fna";
    exit 1;
fi

Threads="$1";
MinPen="$2"
ClustInfo="$3";
Tree="$4";
POInfo="$5";
ClustAssign="$6"
shift 6;


AssignDir="/home/zac/scripts/Hendrik/GutMicrobiome/Round2";
scriptDir="/home/zac/scripts/Hendrik/COVID";
perl "$AssignDir/AssignGene.pl" -vt "$Threads" "$ClustInfo" "$Tree" > "$ClustAssign"
if [ "$?" -ne 0 ]; then >&2 echo "Error in Assigning clusters to txid"; exit 1; fi
perl "$scriptDir/CreatePseudoOrganisms.pl" -p "$MinPen" -o "$POInfo" "$ClustAssign" "$@";
if [ "$?" -ne 0 ]; then >&2 echo "Error writing clusters to pseudo-orgs"; exit 1; fi
