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
  * MAFFT version 7.0+
  * MUSCLE version 3.8.425+
  * CLUSTAL OMEGA version 1.2.1+
* An LCR masker (can be disabled), the following are supported by default:
  * dustmasker version 1.0.0+
  * sdust version 1.0.0+
* BLASTn version 2.9.1+

## Outputs
* A fasta formatted set of probe sequences
* A tab delimited file describing the source and position of the probe sequences
* A tab delimited file detailing the composition of gene clusters
* A fasta formatted set of cluster sequences
* A file containing a newick formatted tree describing the relationship between taxa


## Installation

1. Download the repository
2. Navigate into the master directory
3. Alter the ALN and LCR variables in the Makefile if you wish to use non-default dependencies
4. Run `make`

## Quick Guide

Given that you have a tab delimited file with columns of paths-to-annotated-genome and taxon-id
The pipeline can be run with 
`bin/HUBDesign GenomeInfo.tab`

Further options for the pipeline are detailed with the --help option 
`bin/HUBDesign --help`

## Optional inputs

* `--guide-tree` 
A tree which can be used to guide how the nested relationships between the input genomes. 
By default the hierarchy is determined by the co-occurrence of clusters between taxa.

* `--blast-db` 
Blast databases of blacklist sequences can be provided with the option either with multiple
uses of the flag or with a comma separated list

* `--probe-count` 
The maximum size of the output probe set. By default the maximum number of probes will be a multiple of
the number of input genomes

* `--probe-length`  
The length of probes to design


### Common Parameters

* `--r2t-divergence` 
The maximum amount of divergence within a cluster. As more divergent baits are less effective, selecting
an appropriate value here allows for breadth of coverage while maintaining specificity

* `--penetrance`
The minimum proportion of a node's descendants which have a particular gene for that gene to be included
in the pseudo-genome for that node. High values potentially allow for non-identifying horizontally transferred elements, while low values eliminate many regions from consideration for cross-reactivity reducing specificity of the final probe set.

* `--tiling-density`
The minimum tiling density to aim for during probe selection. Higher tiling density results in more
effective capture of targets, but setting the target too high can make it impossible to balance the
number of baits across organisms. This can cause unexpectedly-small probe sets


### Configuration

Altering the HUBDesign.cfg, or providing a different config file will change the default parameters of the pipeline for ease of use with multiple profiles

The following can be used to generate config files for editing
`bin/Config.pl`

## Tutorial

A directory of test files is provided in the repository:
* Genome files were created by running PROKKA on the raw genome files with the kingdom set to viruses,
  names of the files are the taxon-ids of the genomes
* The GenomeInfo file was created using awk to print out the path and the basename of each file
* The Guide.tree file was generated from the lineages provided in NCBI's taxonomy for 56 refseq
  coronavirus genomes

The following will generate a directory with all final and most intermediate files HUBDesign produces:  
`bin/HUBDesign --guide-tree test/Guide.tree --tiling-density 1 --probe-count 2500 --output-dir test --verbose --keep test/GenomeInfo.tab 2>test.log` 

The output of this command can be compared to the output provided in the test directory

## Advanced pipeline use

The individual phases of the pipeline can be run separately, which is useful if you would like to do
additional processing at any given step. (Ex. Filtering against NCBI's non-redundant nucleotide database for non-viral hits) 
As long as formats are preserved between steps this is relatively painless


### Clustering
`bin/Clustering.pl [-options] GFF_Files ... > RepSeq.fna`

### Assignment
`bin/Assignment.pl [-options] ClusterInfo ... > Assignment.tsv`

### Identification
`Identification.pl [-options] Assignment.tsv RepSeq.fna > CandidateBaitRegions.tab`

### Filtering
`Filtering.pl [-options] Candidates.tab [BlastDB_1 ...] > FilteredCandidates.tab`

### Selection
`Selection.pl [-options] -n maxBaits Tree Candidates > TiledBaits.fna`

