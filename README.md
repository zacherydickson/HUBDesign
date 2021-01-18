# HUBDesign
Hierarchical Unique Bait Design for simultaneous and specific capture of known and novel targets

## Description
Given a set of annotated genomes and a species tree creates a probe set which will allow unique capture of as many nodes in the tree as possible,
attempting to balance the number of baits per genome.

## Inputs
* Annotated Genomes - PROKKA is a suggested annotation tool
* A tree of the input genomes, this will be used to guide clustering

## Requirements
This version requires
* Perl version 5.10.0+
* MAFFT version 7.0+
* distmat from the EMBOSS package
* neighbor from the PHYLIP package
* sdust
* blastn
* awk

## Outputs
A fasta formated set of probe sequences
A tab delimited file describing the source and position of the the probe sequences

## Notes
A User friendly and comprehensive version of the program is a work in progress

Clustering can be performed with a single command:
    RunClustering.pl [-options] GFFFiles ... > RepSequences.fna
