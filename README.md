5/7/2015

This set of scripts is freely available under the GNU license. However, there is no such thing as free lunch :)
These scripts are made available without any guarantees. However, I will try to help wherever I can. Do not hesitate to shoot me an email at koenraad{dot}vandoorslaer{at}gmail{dot}com.


This set of scripts is used to annotate papillomavirus genomes for the PaVE database (http://pave.niaid.nih.gov).
These scripts have been tested on python2.7 on Mac Os X. However, I assume that any python 2.x will work.

Briefly, the program does the following steps (for more detail, see below)

1) download GenBank file
2) translate the Sequence in three forward open reading frames
3) Use a custom Blast database to annotate viral proteins
4) Attempt to annotate E1^E4 and E8^E2 spliced transcripts
5) look for consensus E2 binding sites in the genomes
6) Identify the viral URR
7) perform some inital QC steps

Requirements
I apologize for the many dependencies. Contact me if installation is not going as expected.
1) Python
2) Biopython (http://biopython.org/DIST/docs/install/Installation.html)
3) Mechanize (http://wwwsearch.sourceforge.net/mechanize/download.html)
4)ete2 (https://pythonhosted.org/ete2/install/)
5)stand-alone BLAST (http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/)
      do not forget to initialize the BLAST database (see below)

Detailed steps

3) Use a custom Blast database to annotate viral proteins
    The file "DB.blast.txt" contains the data to create a custom blast database. It contains the example protein          sequences for each of the viral ORFs. For E5, I chose to use consensus sequences to improve the blast annotation      performance. Furthermore, the database includes many random sequences to improve the statistics of the BLAST          search. Take the following steps to initialize the DB.
        1) cd into the "Blast" folder  
        2) type (without the quoutes) "makeblastdb DB.blast.txt"
    
    The script will translate all ORFs in the forward frames and compare these putative proteins to the Blast             database. If a hit is found, this ORF will be annotated based on the Blast result.
    Annotation of L1 tries to take into account that L1 is usually translated from a spliced mRNA. The regular            expression "(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG" will ensure that the correct methionine is used to start the L1 ORF.

4) Attempt to annotate E1^E4 and E8^E2 spliced transcripts
	The annotation process leverages a published server to predict aplice donors and annotators (Wang M. and Marín A. 2006. Characterization and Prediction of Alternative Splice Sites. Gene 366: 219-227.). The annotion script uses the default settings as described in the original paper. The result table is parsed based on several assumptions. These assumptions are the result of extensive comparative genomics between well studies papillomaviruses.
	E1 splice donor assumptions:
		1) The SD is located within the E1 ORF
		2) The canonical "AG" dinucleotide has to be used
		3) When spliced into the E4/E8 splice acceptor, the resulting protein has to be in frame with the 			   annotated E4 ORF
	E8 splice donor assumptions:
		1) The E8 ORF is located in the "+1 frame" of E1
		2) The canonical "AG" dinucleotide has to be used
		3) When spliced into the E4/E8 splice acceptor, the 
		resulting protein has to be in frame with the annotated E2 ORF
	E4/E8 splice acceptor assumptions:
		1) E4 is loacted in the "+1 frame" of E2
		2) The canonical "AG" dinucleotide has to be used

Please Note: while the process works fairly well, the produced spliced sequences are manually compared to sequences currently in PaVE. When needed, manual curration is performed.

5) look for consensus E2 binding sites in the genomes
The relaxed consensus (ACCnnnnnnGGT) is used to search the viral genome. Importantly, the script takes into account that the linear sequence, in fact represents a circular genome.

6) Identify the viral URR
Traditionally, the URR is located between the end of L1 and the beginning of the next annotated ORF

7) perform some inital QC steps
	1) test that the coordinates of the different ORFs are correct.
		E.g. E1^E4 should end at the same coordinate as E4
		E.g. all ORFs should be located within the viral genome
	2) Check whether expected ORFs are annotated
	Based on comparative genomics it is expected that viral genomes contain certain ORFs. These requirements vary 	depending on the phylogenetic classification of these viruses. E.g. members of the Xipapillomavirus genus do not 		contain an E6 ORF. For this step, the L1 sequence is extracted and compared to all L1 sequences currently 	in PaVE. This comparison is used to classify the test genome into a viral genus and/or species. Based on this 		certain recommendations are made. This step also ensures that no ORF is duplicated. E.g. A virus cannot 		contain multiple L1 ORFs. The presence of such a duplication usually indicates a (sequencing) error in the 		genome. 
The results of these QC steps is stored in the log.txt file





