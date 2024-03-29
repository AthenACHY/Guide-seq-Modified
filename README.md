# Guide-seq-Modified
We used modified Guide-seq codes for off-target detections of SaCas9 variants. The original code are posted in [guideseq: The GUIDE-Seq Analysis Package](https://github.com/aryeelab/guideseq).

## Details
identifyOfftargetSites_1.py has its mismatches evaluation between the sgRNA template and the target sequence changed to acommodate the flexible SaCas9 PAM site "NNGRRT".


```
python /home/athena/bin/guideseq/guideseq/identifyOfftargetSites_SA.py \
--samfile ${outfile}/aligned/${readID[i]}.sam \
--ref /data/genome_reference/Homo_sapiens_assembly19.fasta \
--outfile ${outfile}/identified/${case[i]}_identifiedOfftargets.txt \
--target ${Guide[i]} --mismatch 7
```


visulization_SA_pdfout.py is changed to allow the display of the nucleotide "R" in the off-target panel and generate a pdf file for the figure in addition to the .svg format.  


```
python /home/athena/bin/guideseq/guideseq/visualization_SA_pdfout.py \
${outfile}/identified/${case[i]}_identifiedOfftargets.txt \
${outfile}/visualise/${case[i]}_identifiedOfftargets $title
```
## Additional codes
2018_06_20_Amp_analysis.py is a in-house script used for analysing editing efficiency of SaCas9 from the amplicon ssequencing data. 

2019_01_10_characterize_indel_type.py characterize per base insertion/deletion consequential of Cas9 cleaveage and repair of the VEGFA_8 target site.
