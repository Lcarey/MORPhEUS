# MORPhEUS
RNA polymerase errors cause splicing defects and can be regulated by differential expression of RNA polymerase subunits

Errors during transcription may play an important role in determining cellular phenotypes: the RNA polymerase error rate is >4 orders of magnitude higher than that of DNA polymerase and errors are amplified >1000-fold due to translation. However, current methods to measure RNA polymerase fidelity are low-throughout, technically challenging, and organism specific. Here I show that changes in RNA polymerase fidelity can be measured using standard RNA sequencing protocols. I find that RNA polymerase is error-prone, and these errors can result in splicing defects. Furthermore, I find that differential expression of RNA polymerase subunits causes changes in RNA polymerase fidelity, and that coding sequences may have evolved to minimize the effect of these errors. These results suggest that errors caused by RNA polymerase may be a major source of stochastic variability at the level of single cells.

Lucas B. Carey. eLife 2015. 

-------------------

mpileupCountMismatches.pl takes as input the output of samtools mpileup, and produces a tab file with the number of errors and number of ref reads at each position in the genome. The overall error rate is the sum of the last column divided by the sum of the second to last column. 

The Makefile provides an example Makefile to download sequencing reads from NCBI, run bwa, and run MORPhEUS


