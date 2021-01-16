#!/bin/bash

Usage="$(basename $0) [-options] RegionInfo > FilteredRegionInfo"

windowLen=64;
threshold=20;
while getopts w:t: arg; do
	case $arg in
		#OPTIONS
		w)
			windowLen=$OPTARG
			;;
		t)
			threshold=$OPTARG
			;;
		#FLAGS
		#HELP
		?)
			>&2 echo -e " ===Description===\n" \
				"===Usage===\n" \
                "\t$Usage\n" \
				"===Parameters===\n" \
                "\tRegionInfo\tA tab delimited file with 4 columns: RegID, target, position, and sequence\n" \
				"===Options===\n" \
				"\t-w INT\tWindow Length (Default: 64)\n" \
				"\t-t INT\tComplexity Threshold(Default: 20)\n" \
				"===Flags===\n" \
				"\t-?\t Display this message and exit\n;"
			exit
		;;
	esac
done
shift $((OPTIND-1));


if [ "$#" -lt 1 ]; then
    >&2 echo -e "Usage: $Usage\n\tUse -? for more info";
    exit 1;
fi;

RegionInfoFile=$1;

sdust <(awk -F '\\t' '{print ">"$1; print $4}' "$RegionInfoFile")  | 
    awk -v minRegLen="75" '
        (ARGIND == 1){
            EntryCount[$1]++;
            LCRPos[$1"s"EntryCount[$1]] = $2;
            LCRPos[$1"e"EntryCount[$1]] = $3
            next;
        }
        {
            if(!EntryCount[$1]){print $0;next}
            init=0;
            for(i=1;i<=EntryCount[$1];i++){
                if(LCRPos[$1"s"i] >= init+minRegLen){
                    print $1"_"i"\t"$2"\t"$3+init"\t"substr($4,init+1,LCRPos[$1"s"i]-init);
                }
                init = LCRPos[$1"e"i]+1;
            }
            if(length($4) - init >= minRegLen){
                print $1"_"EntryCount[$1]+1"\t"$2"\t"$3+init"\t"substr($4,init+1);
            }
        }
        ' /dev/stdin "$RegionInfoFile"
