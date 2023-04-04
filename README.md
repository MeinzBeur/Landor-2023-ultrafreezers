# Ten degrees greener: DNA, RNA, and prokaryote community sample stability at different ultra-low temperature storage conditions

Lotta A. I. Landor, Thomas M. J. Stevenson, Kyle M. J. Mayers, Mitchell S. Fleming, Sven Le Moine Bauer, Hannah Babel, Stefan Thiele

This repository contains the scripts and files needed to reproduce the data analysis presented in the aforementioned article. Note that the authors are in no case Unix/R professionals, and the code can certainly be written in a more idiomatic way. Do not hesitate to reach out for further help.

The following links will bring you to:
- The compositional analysis: [Figure 3A and Permanova test](CoDA-analysis.md) 
- The dispersion analysis: [Supplementary data 4](SuppMaterial_Dispersion.md)

The following files are also given:
- [Metadata.csv](Metadata.csv): The context information to the community samples.
- [Otutab.sorted.tsv](Otutab.sorted.tsv): The OTU table produced by the sequence processing pipeline.
- [assignments.csv](assignments.csv): The taxonomy information assigned by CREST4, prior to decontamination.
- [Otus_curated.fasta](Otus_curated.fasta): Centroids of the OTUs. Not used in any script presented here.

We do not describe here the making of figure 1 and 2, as they do not involve any analysis, just plotting. For the barplots presented in Supplementary material 3, we refer to the script described [there](https://github.com/MeinzBeur/LeMoineBauer-2022-Milos/blob/main/Barplots.md).

The pipeline describing the processing of the sequences prior to statistical analysis can be found [there](https://github.com/MeinzBeur/LeMoineBauer-2022-Milos/blob/main/Pipeline%20explanations.md). The fastq files can be recovered from the European Nucleotide Archive under the project number PRJEB60555. 
