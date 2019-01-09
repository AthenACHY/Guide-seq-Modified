# Guide-seq-Modified
We used modified Guide-seq codes for off-target detections of SaCas9 variants. The original code are posted in [guideseq: The GUIDE-Seq Analysis Package](https://github.com/aryeelab/guideseq).

## Details
identifyOfftargetSites_1.py has its mismatches evaluation between the sgRNA template and the target sequence changed to acommodate the flexible SaCas9 PAM site "NNGRRT".

visulization_SA_pdfout.py is changed to allow the display of the nucleotide "R" in the off-target panel and generate a pdf file for the figure in addition to the .svg format.  

2018_06_20_Amp_analysis.py is a in-house script used for analysing editing efficiency of SaCas9 from the amplicon ssequencing data. 

