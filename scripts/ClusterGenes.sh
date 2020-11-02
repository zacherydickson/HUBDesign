#!/bin/bash

if [ "$#" -lt 6 ]; then
    >&2 echo -e "Usage: $(basename $0) Threads GeneticCode Root2TipDist ClustInfoOut GeneInfo.tab " \
    "SequenceFiles > ConsensusClusters.fna\n" \
        "\tNote: Sequence Files should be provided in the same order as for GeneInfo"
    exit 1;
fi


Threads="$1";
GenCode="$2";
R2TDist="$3";
ClustOut="$4";
shift 4;

WorkingDir="/scratch/zac/clustgene_$(rndstr 8)";
export TMPDIR="$WorkingDir";
>&2 echo "WorkingDir: $WorkingDir";
mkdir -p "$WorkingDir";

##Step One get gene files
scriptDir="/home/zac/scripts/Hendrik/COVID";
perl "$scriptDir/ExtractGenes.pl" "$WorkingDir" "$@";
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in Extract Genes";
    exit 1;
fi


##Deal with any single sequences
ls "$WorkingDir/"*".fna" | 
    parallel -j "$Threads" 'if [ $(grep -c "^>" {}) -lt 2 ]; then sed -i "s/>.\+/>Cluster_0/" {};mv {} {.}.consensus; echo -e "1\t0" > {.}.clust; fi'
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in Handling Lone Sequences";
    exit 1;
fi

#Get an MSA for each gene
#Run Clustering on each Gene

for f in "$WorkingDir/"*".fna"; do
    sample=$(basename "$f" ".fna");
    if [ ! -f "$WorkingDir/$sample.msa" ]; then 
        "$scriptDir/AlignGenesByProt.sh" "$Threads" "$GenCode" "$f" > "$WorkingDir/$sample.msa";
        if [ "$?" -ne 0 ]; then
            >&2 echo "ERROR Alignment by Prot for $f";
            exit 1;
        fi
    fi
done

export ClusterDir="/home/zac/scripts/Hendrik/GutMicrobiome/Round2"
export R2TDist="$R2TDist";
ls "$WorkingDir/"*".msa" | parallel -j "$Threads" '"$ClusterDir/ClusterAlignment.sh" $R2TDist {.}.clust {}'
if [ "$?" -ne 0 ]; then
    >&2 echo "ERROR in Clustering Alignments";
    exit 1;
fi

##Generate ClusterInfo
perl "$scriptDir/CreateClusterInfo.pl" "$1" "$WorkingDir/"*".clust" > "$ClustOut";

rm -rf $WorkingDir
