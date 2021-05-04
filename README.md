# HUBDesign
Hierarchical Unique Bait Design for simultaneous and specific capture of known and novel targets

## Description
Given a set of annotated genomes and a species tree creates a probe set which will allow unique capture of as many nodes in the tree as possible,
attempting to balance the number of baits per genome.

## Inputs
* Annotated Genomes - PROKKA is a suggested annotation tool
* A tree of the input genomes, this will be used to guide clustering
* A set of blast databases of blacklist sequences

## Requirements
This version requires
* Perl version 5.10.0+
* BioPerl 1.7.10+
* A multiple sequence aligner, the following are supported by default:
.. * MAFFT version 7.0+
.. * MUSCLE version 3.8.425+
.. * CLUSTAL OMEGA version 1.2.1+
* An LCR masker (can be disabled), the following are supported by default:
.. * dustmasker version 1.0.0+
.. * sdust version 1.0.0+
* BLASTn version 2.9.1+

## Outputs
* A fasta formated set of probe sequences
* A tab delimited file describing the source and position of the the probe sequences

* A tab delimited file detailing the composition of gene clusters
* A fasta formated set of cluster sequences
* A file containing a newick formatted tree describing the relationship between taxa


## Installation

1. Download the repository
2. Navigate into the master directory
3. Compile the SA_BOND executable
.. 1. Navigate into the SA_BOND directory
.. 2. `g++ -fopenmp src/sabond.cpp -o sabond`
4. Set up Default Configs
.. 1. Return to root of master directory
.. 2. `bin/Configure -a multiple-sequence-aligner-of-choice -l lcr-masker-of-choice > HUBDesign.cfg`

## Quick Running the pipeline

Given that you have a tab delimited file with columns of paths-to-annotated-genome and taxon-id
The pipeline can be run with 
`bin/HUBDesign GenomeInfo.tab`

A tree informaing the relationships between genomes can be provided with the --guide-tree option

Blast databases of blacklist sequences can be provided with the --blast-db option either with multiple
uses of the flag or with a comma separated list

Further options for the pipeline are detailed with the --help option

Altering the HUBDesign.cfg, or providing a different config file will change the default parameters of the pipeline for ease of use with multiple profiles

The foollowing can be used to generate config files for editing
`bin/Config.pl`

##Tutorial

A directory of test files is provided in the repository

Running `bin/HUBDesign --verbose --keep test/GenomeInfo.tab 2>test.log` 

Will generate a directory with all final and intermediate files HUBDesign produces


## Advanced pipeline Use

The individaul phases of the pipeline can be run separately, which is useful if you would like to do
additional processing at any given step. (Ex. Filtering against NCBI's non-redundant nucleotide database for non-viral hits)
As long as formats are preserved between steps this is relatively painless


###Clustering
`bin/Clustering.pl [-options] GFF_Files ... > RepSeq.fna`

###Assignment
`bin/Assignment.pl [-options] ClusterInfo ... > Assignment.tsv`

###Identification
`Identification.pl [-options] Assignment.tsv RepSeq.fna > CandidateBaitRegions.tab`

###Filtering
`Filtering.pl [-options] Candidates.tab [BlastDB_1 ...] > FilteredCandidates.tab`

###Selection
`Selection.pl [-options] -n maxBaits Tree Candidates > TiledBaits.fna`

