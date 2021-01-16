
============================================================
How to run SA_BOND (our program for designing oligonucleotides)
============================================================

	./SA_BOND <inputSequences> <outputOligos> [-length]
[-seqSim] [-maxMatch] [-maxGC] [-minGC] [-dimerLen] [-dimerSim] [-hairpinLen]
[-minhpLoop] [-maxhpLoop] [-rangeTm] [-minTm] [-maxTm] [-oligCon] [-saltCon]
[-secStr] [-outMode] [-simCheck] [-paddingLen]
[-tilingLen] [-maxThread]

Required parameters:

<inputSequencess> - the input file containing the gene sequences in FASTA format
<outputOligos> - the output file containing the designed oligos and additional information (gene number, oligo sequence, length, melting temperature, distance from 5' and 3' ends)


Optional parameters:

[-length]:     The length of oligo (Default value is 50)
[-seqSim]:     Sequence identity percentage (Default value is 75)
[-maxMatch]:   Maximum consecutive match (Default value is 15)
[-maxGC]:      Maximum GC-content percentage (Default value is 70)
[-minGC]:      Minimum GC-content percentage (Default value is 30)
[-dimerLen]:   The length of dimer window (Default value is 15)
[-dimerSim]:   Dimer stringency percentage (Default value is 85)
[-hairpinLen]: The length of hairpin stem (Default value is 6)
[-minhpLoop]:  Minimum length of hairpin loop (Default value is 1)
[-maxhpLoop]:  Maximum length of hairpin loop (Default value is 3)
[-minTm]:      Minimum melting temperature (Disabled in default mode) 
[-maxTm]:      Maximum melting temperature (Disabled in default mode)
[-rangeTm]:    The length of melting temperature interval (Default value is 10)
[-oligCon]:    Oligo concentration in nM (Default value is 1000)
[-saltCon]:    Salt concentration in mM (Default value is 75)
[-secStr]:     Enabling secondary structure checking
[-maxOligs]:   The number of oligos to return per input sequence (Default of 1)
                    A value of -1 returns all valid oligos
[-outMode]:    The oligo output mode 0: (default) print each oligo and its info on its own line
                                mode 1: collapse continuous oligos together onto one line (maxOligs must be -1)
[-simCheck]:   Whether or not to perform a fast homology check (useful reduce search space on an exhaustive search)
                    Note: This similarity check is not strand aware
[-paddingLen]: The amount of flanking DNA for each oligo on both the 5 and
3 prime ends (Default value is 0).
                Note: This is a maximum value (oligos near the edges of genes
may have less padding)
[-tilingLen]:  The length of baits to be tiled within the olig (Default value
is the same as -length)
                Note: Output sequences will be regions of up to "-length" bp
made up of up to "-length" - "-tilingLen" + 1 unique "-tilingLen" oligos
[-maxThread]: The maximum number of threads to use (Default value is system max)


Example
=======

For the input dataset "ecoli.fsa" the commands for running BOND with default parameters on linux, mac, and windows are:

1. For Linux:
	./BOND_linux ecoli.fsa ecoli.bond

2. For Mac:
	./BOND_mac ecoli.fsa ecoli.bond

3. For Windows (64-bit):
	BOND_windows ecoli.fsa ecoli.bond


If the user would like to design oligos with these parameters:

- length: 60
- sequence identity percentage: 80
- maximum consecutive match: 16
- enabling secondary structure checking

then the linux command is:

	./BOND_linux ecoli.fsa ecoli.bond -length 60 -seqSim 80 -maxMatch 16 -secStr 



=============================================================================================================
How to run EVAL (the program for evaluating the specificity); PRE_EVAL prepares the oligo file for evaluation
=============================================================================================================

	./PRE_EVAL_<platform> <BONDOligos> <plainOligos>
	./EVAL_<platform> <inputSequences> <plainOligos> <outputReport>

Required parameters:

<BONDOligos> - the file containing the oligos designed by BOND
<plainOligos> - a file containing only the oligos, one per line, and no additional information
<inputSequencess> - the input file containing the gene sequences in FASTA format
<outputReport> - the output file containing the evaluation report (bad oligos, good oligos, and statistics)


Example
=======

To evaluate the oligos designed by BOND "ecoli.bond" for the dataset "ecoli.fsa" the linux commands are:

        ./PRE_EVAL_linux ecoli.bond ecoli.plain
        ./EVAL_linux ecoli.fsa ecoli.plain ecoli.report



=====================================================================================================
How to run MAX (the program for estimating the maximum number of oligos, needed to evaluate coverage)
=====================================================================================================

	./MAX_<platform> <inputSequences> <outputReport> 

Required parameters:

<inputSequencess> - the input file containing the gene sequences in FASTA format
<outputReport> - the output file containing the estimated maximum number of oligos


Example
=======

To estimate the maximum number of oligos for the dataset "ecoli.fsa" the linux command is:

        ./MAX_linux ecoli.fsa ecoli.max




=================================================
For any questions, please contact 
Lucian Ilie 
e-mail: ilie@csd.uwo.ca
=================================================

