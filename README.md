# nucleolus-centromere-tsa-data-analysis

## This is the code used to analyze the data presented in Pradeep Kumar et al., 2023, titled "Nucleolus and centromere TSA-Seq reveals variable localization of heterochromatin in different cell types". Here is the link to the preprint: https://doi.org/10.1101/2023.10.29.564613

This folder has three sub folders:

Folder 1: histone_mark_overlap <br />
The R code in this folder was used for calculating overlap of histone marks, repli-seq and RNA-seq with different BED file for example LAD subset 1. This code is written for one cell line but this code was edited to change cell line name and BED file for analysis. The output of this code can be processed in matlab or R to create bar graphs. This code was used for analysis shown in Fig 5, Fig 6, Fig S5 and Fig S7.  

Folder 2: Mean_fish_distance_calculation <br />
This Matlab code was used to calculate the mean distance of FISH probes from Su, J.H., et al., "Genome-Scale Imaging of the 3D Organization and Transcriptional Activity of Chromatin". Cell, 2020. 182(6): p. 1641-1659 e26.
The processed excel sheet (fish.xlsx) of the FISH data from Su et al., is also provided. The output of this code can be processed in matlab or R to create scatter plots. The mean data from this FISH analysis was used in Fig 3e, 3h, 3k and Fig 7g.  

Folder 3: Type1_2_peaks_vs_TSA_Seq <br />
This R code was used to overlap Type 1 and Type 2 peaks with TSA-seq data. The output of this code can be processed in matlab or R to create scatter plots. This code was used in Fig 7b and 8b.


