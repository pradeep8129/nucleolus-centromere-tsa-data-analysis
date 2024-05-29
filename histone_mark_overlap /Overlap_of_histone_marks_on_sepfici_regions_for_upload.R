rm(list=ls())

# Set working directory that contains all the files 
setwd("~path to the folder")   

# Load relevant libraries

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(zoo)
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(data.table)
library("readxl")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(regioneR)
library(readr)

####### importing LAD calls for HCT116 and HFFc6 cells #######

h1_damid_calls <- import("H1_LMNB1-25kb-combined_AD.bed")   # H1 LAD calls
hct_damid_calls <- import("HCT116_LMNB1-25kb-combined_AD.bed")  # HCT116 LAD calls
hg38_100kb <- import("hg38.100kb.bed")  # (hg38 bins)


###### Identifying bins for LADs and iLADs in HCT116 #####

hg38_100kb$status = "ilAD"
hct_lad_subjects <-queryHits(findOverlaps(query = hg38_100kb, subject = hct_damid_calls))
hg38_100kb [hct_lad_subjects]$status = "LAD"

hct_lad <- granges(hg38_100kb[hg38_100kb$status == "LAD"])
hct_ilad <- granges(hg38_100kb[hg38_100kb$status == "ilAD"])

######### Identifying bins for LADs and iLADs in H1 ######

hg38_100kb_copy <- granges(hg38_100kb)
hg38_100kb_copy$status = "ilAD"
h1_lad_subjects <-queryHits(findOverlaps(query = hg38_100kb_copy, subject = h1_damid_calls))
hg38_100kb_copy[h1_lad_subjects]$status = "LAD"

h1_lad <- granges(hg38_100kb_copy[hg38_100kb_copy$status == "LAD"])
h1_ilad <- granges(hg38_100kb_copy[hg38_100kb_copy$status == "ilAD"])

############### Importing ChIP seq data of H1 and HCT116 cells #####

h1.chip<-read_tsv("H1_chip_anno_25kb_adjusted_no_na.txt")         # H1 ChIP seq data
hct.chip <-read_tsv("HCT116_chip_anno_25kb_adjusted_no_na.txt")   # HCT116 ChIP seq data

############### Importing Specific bed file for example LAD subset 1 #####

fig_5a_bed_file <- import("LAD Subset1.bedfig_5a_new_bed")   # Bed file
high_mki_h1_low_mki_hct116 <- granges(fig_5a_bed_file)

######## Making granges from tsv files of ChIP seq for H1 and HCT116 #####

h1.chip.gr <- makeGRangesFromDataFrame(h1.chip,keep.extra.columns = T, ignore.strand=FALSE,
                                       seqinfo=NULL,
                                       seqnames.field=c("seqnames", "seqname",
                                                        "chromosome", "chrom",
                                                        "chr", "chromosome_name",
                                                        "seqid"),
                                       start.field="start",
                                       end.field=c("end", "stop"),
                                       strand.field="strand",
                                       starts.in.df.are.0based=FALSE)


start(h1.chip.gr)<-start(h1.chip.gr)+1


hct.chip.gr <- makeGRangesFromDataFrame(hct.chip,keep.extra.columns = T, ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames", "seqname",
                                                         "chromosome", "chrom",
                                                         "chr", "chromosome_name",
                                                         "seqid"),
                                        start.field="start",
                                        end.field=c("end", "stop"),
                                        strand.field="strand",
                                        starts.in.df.are.0based=FALSE)


start(hct.chip.gr)<-start(hct.chip.gr)+1

##### Finding overlaps between histone marks and specific bed regions for example Subset 1 ######
hct_chip_high_mki_h1_low_mki_hct116 <- hct.chip.gr[unique(queryHits(findOverlaps(hct.chip.gr,high_mki_h1_low_mki_hct116)))]
h1_chip_high_mki_k562_low_mki_hct116 <- h1.chip.gr[unique(queryHits(findOverlaps(h1.chip.gr,high_mki_h1_low_mki_hct116)))]

#### exporting csv files #### 

write.table( x = data.frame(hct_chip_high_mki_h1_low_mki_hct116), file = "overlap_hct_chip_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_chip_high_mki_k562_low_mki_hct116), file = "overlap_hff_chip_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

##### Finding overlaps between Repli-Seq and specific bed regions for example Subset 1 ######

hct_repli <- import("HCT116.smooth.bw")     # HCT repli seq
h1_repli <- import("HFF.smooth.bw")         # H1 repli seq

hct_repli_high_mki_h1_low_mki_hct116 <- hct_repli[unique(queryHits(findOverlaps(hct_repli,high_mki_h1_low_mki_hct116)))]
h1_repli_high_mki_h1_low_mki_hct116 <- h1_repli[unique(queryHits(findOverlaps(h1_repli,high_mki_h1_low_mki_hct116)))]

#### exporting csv files #### 

write.table( x = data.frame(hct_repli_high_mki_h1_low_mki_hct116 ), file = "overlap_subset__hct_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_repli_high_mki_h1_low_mki_hct116 ), file = "overlap_subset_hff_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

####### Calculating the histone marks for LADs and iLADS in HCT116 and H1 #######
hct_lads_histone_high_mki_h1_low_mki_hct116 <- hct.chip.gr[unique(queryHits(findOverlaps(hct.chip.gr,hct_lad)))]
hct_interlads_histone_high_mki_h1_low_mki_hct116 <- hct.chip.gr[unique(queryHits(findOverlaps(hct.chip.gr,hct_ilad)))]

h1_lads_histone_high_mki_h1_low_mki_hct116 <- h1.chip.gr[unique(queryHits(findOverlaps(h1.chip.gr,h1_lad)))]
h1_interlads_histone_high_mki_h1_low_mki_hct116 <- h1.chip.gr[unique(queryHits(findOverlaps(h1.chip.gr,h1_ilad)))]

####### Calculating the replication timing for LADS and iLADS in HCT116 and H1 #######
hct_lads_repli_high_mki_h1_low_mki_hct116 <- hct_repli[unique(queryHits(findOverlaps(hct_repli,hct_lad)))]
hct_interlads_repli_high_mki_h1_low_mki_hct116 <- hct_repli[unique(queryHits(findOverlaps(hct_repli,hct_ilad)))]

h1_lads_repli_high_mki_h1_low_mki_hct116 <- h1_repli[unique(queryHits(findOverlaps(h1_repli,h1_lad)))]
h1_interlads_repli_high_mki_h1_low_mki_hct116 <- h1_repli[unique(queryHits(findOverlaps(h1_repli,h1_ilad)))]

#### exporting csv files #######

write.table( x = data.frame(hct_lads_histone_high_mki_h1_low_mki_hct116), file = "overlap_hct_lads_histone_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(hct_interlads_histone_high_mki_h1_low_mki_hct116), file = "overlap_hct_interlads_histone_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

write.table( x = data.frame(h1_lads_histone_high_mki_h1_low_mki_hct116), file = "overlap_hff_lads_histone_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_interlads_histone_high_mki_h1_low_mki_hct116), file = "overlap_hff_interlads_histone_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )


write.table( x = data.frame(hct_lads_repli_high_mki_h1_low_mki_hct116), file = "overlap_hct_lads_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(hct_interlads_repli_high_mki_h1_low_mki_hct116), file = "overlap_hct_interlads_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

write.table( x = data.frame(h1_lads_repli_high_mki_h1_low_mki_hct116), file = "overlap_hff_lads_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_interlads_repli_high_mki_h1_low_mki_hct116), file = "overlap_hff_interlads_repli_high_mki_h1_low_mki_hct116.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

############### Importing RNA Seq data for HCT116 and H1  #####
h1.chip_rna_Seq<-read_excel("h1_expression.xlsx")
hct.chip_rna_seq<- read_excel("hct_expression.xlsx")

######## Creating granges from tsv files #####

h1.chip_rna_Seq_gr <- makeGRangesFromDataFrame(h1.chip_rna_Seq,keep.extra.columns = T, ignore.strand=FALSE,
                                               seqinfo=NULL,
                                               seqnames.field=c("seqnames", "seqname",
                                                                "chromosome", "chrom",
                                                                "chr", "chromosome_name",
                                                                "seqid"),
                                               start.field="start",
                                               end.field=c("end", "stop"),
                                               strand.field="strand",
                                               starts.in.df.are.0based=FALSE)


hct.chip_rna_Seq_gr <- makeGRangesFromDataFrame(hct.chip_rna_seq,keep.extra.columns = T, ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c("seqnames", "seqname",
                                                                 "chromosome", "chrom",
                                                                 "chr", "chromosome_name",
                                                                 "seqid"),
                                                start.field="start",
                                                end.field=c("end", "stop"),
                                                strand.field="strand",
                                                starts.in.df.are.0based=FALSE)

##### inding overlaps between RNA-seq and specific bed regions for example Subset 1 ######

h1_fig5a_rna_seq <- h1.chip_rna_Seq_gr[unique(queryHits(findOverlaps(h1.chip_rna_Seq_gr,high_mki_h1_low_mki_hct116)))]
hct_fig5a_rna_seq <-hct.chip_rna_Seq_gr[unique(queryHits(findOverlaps(hct.chip_rna_Seq_gr,high_mki_h1_low_mki_hct116)))]

#### exporting csv files ####

write.table( x = data.frame(hct_fig5a_rna_seq), file = "subset_hct_expression.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_fig5a_rna_seq), file = "subset_hff_expression.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

####### now calculating the histone marks for lads #######
hct_lads_fig5<- hct.chip_rna_Seq_gr[unique(queryHits(findOverlaps(hct.chip_rna_Seq_gr,hct_lad)))]
hct_interlads_fig5 <- hct.chip_rna_Seq_gr[unique(queryHits(findOverlaps(hct.chip_rna_Seq_gr,hct_ilad)))]

h1_lads_fig5 <- h1.chip_rna_Seq_gr[unique(queryHits(findOverlaps(h1.chip_rna_Seq_gr,h1_lad)))]
h1_interlads_fig5 <- h1.chip_rna_Seq_gr[unique(queryHits(findOverlaps(h1.chip_rna_Seq_gr,h1_ilad)))]

#### exporting csv files ####

write.table( x = data.frame(hct_lads_fig5), file = "overlap_hct_lads_fig5_rna_seq.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(hct_interlads_fig5), file = "overlap_hct_interlads_fig5_rna_seq.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

write.table( x = data.frame(h1_lads_fig5), file = "overlap_hff_lads_fig5_rna_seq.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(h1_interlads_fig5), file = "overlap_hff_interlads_fig5_rna_seq.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )























