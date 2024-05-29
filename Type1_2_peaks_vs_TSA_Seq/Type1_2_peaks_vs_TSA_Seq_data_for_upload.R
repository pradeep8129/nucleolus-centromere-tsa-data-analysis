# Set working directory that contains all the files 
setwd("~path to the folder")   

# Load relevant libraries

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
#install.packages("zoo")
library(zoo)
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(data.table)
install.packages("readxl")
library("readxl")

############## Importing bw files  ###########

pol1re <- import("HFFc6_Pol1RE_rep1_25kb_conE_hg38_PK_202101.bw")      # bw file of pol1re
mki67ip <- import("HFFc6_MKI67IP_rep1_25kb_conE_hg38_PK_202101.bw")    # bw file of MKI67IP
son <- import("HFFc6_SON_TSA_25kb_hg38_rep1_201806condE.bw")           # bw file of SON

########## Smoothing the tsa-Seq data using  rollapply fuction ##############

#### this is the function I wrote using rollapply to smooth the data ######

bw_moving_window=function(bw,window_size,FUN){
  rol_protein_name<-bw
  #rol_protein_name[rol_protein_name$score==0]$score <-NA #here is converted 0 into NA
  chr_name<-unique(seqnames(rol_protein_name))
  for(i in chr_name){
    rol_protein_name[seqnames(rol_protein_name)==i]$score <-rollapply(rol_protein_name[seqnames(rol_protein_name)==i]$score, width=window_size, FUN=FUN, by= 1, fill = 0, align = "center", partial=TRUE)
    
  }
  return(rol_protein_name)
}

rol_pol1re <- bw_moving_window(pol1re,11,mean)            # here, pol1re is the name of the bw file, 11 is the number of bins used to calculate the mean, and mean is the function applied on the bins
rol_mki67ip <- bw_moving_window(mki67ip,11,mean)
rol_son <- bw_moving_window(son,11,mean)

############## Importing and sorting type 1 and type 2 peaks ###########
type_1_type2 <- import("HFF_TypeI_II_peak.bed")
only_type1 <- granges(type_1_type2[type_1_type2$name == "Type_I_Peaks"])
only_type2 <- granges(type_1_type2[type_1_type2$name == "Type_II_Peaks"])

############## Calculating average TSA-Seq score over the type 1 and type 2 peaks ###########

#### this is the function I wrote to calculate the average score over specific regions provided as bed files ######

bed.avg.score<-function(loc.gr,score.gr,score.col="score"){
  # For every feature in loc.gr this function computes the mean of the score.col in score.gr.
  # the bins you like as GRanges object
  # The GRe=anges that you want to change the binnig
  # name of the column in the score.gr mcols, which will be "score" if you've just imported the file.
  ov <- as.matrix(findOverlaps(loc.gr,score.gr))
  dt=data.table(id=ov[,1],score=values(score.gr)[ov[,2],which(names(values(score.gr))==score.col)])
  dt=dt[,list(av=mean(score)),by=id]
  names(dt)[2]<-paste0(score.col,".mean")
  res=loc.gr[dt$id,]
  mcols(res)=DataFrame(mcols(res),dt[,2])
  res
}


hff_pol1re_type1 <- bed.avg.score(only_type1,rol_pol1re,"score")
hff_mki67ip_type1 <- bed.avg.score(only_type1,rol_mki67ip,"score")
hff_pol1re_type2 <- bed.avg.score(only_type2,rol_pol1re,"score")
hff_mki67ip_type2 <- bed.avg.score(only_type2,rol_mki67ip,"score")
hff_son_type1 <- bed.avg.score(only_type1,rol_son,"score")
hff_son_type2 <- bed.avg.score(only_type2,rol_son,"score")

#### Exporting data as csv files ######

only_type_1_data <- data.frame(Pol1re_Type1_preaks = hff_pol1re_type1$score.mean, MKI67IP_Type1_preaks = hff_mki67ip_type1$score.mean, SON_Type1_peaks = hff_son_type1$score.mean)
write.table( x = data.frame(all_nucleolus_data), file = "HFF_type1.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

only_type_2_data <- data.frame(Pol1re_Type2_preaks = hff_pol1re_type2$score.mean, MKI67IP_Type2_preaks = hff_mki67ip_type2$score.mean, SON_Type2_peaks = hff_son_type2$score.mean)
write.table( x = data.frame(only_type_2_data), file = "HFF_type2.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )


############# Plotting the graph #################

ggplot() +
  geom_point(data = all_nucleolus_data, aes(Pol1re_Type1_preaks,x=SON_Type1_peaks), color = "red")  + # must include argument label "data"
  geom_point(data = only_type_2_data, aes(y=Pol1re_Type2_preaks,x=SON_Type2_peaks)) + geom_smooth(method = "lm")
