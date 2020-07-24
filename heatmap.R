rm(list = ls())
library(readr)
library(plyr)
library(readxl)
library(stringr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### discard pool & iRT & con_
proteinMatrix <- read.table("CVDSSB_proteins_2020426.txt",sep = "\t",header = T,stringsAsFactors = F)

pM <- proteinMatrix[!grepl(pattern ="iRT", proteinMatrix$prot) & !grepl(pattern ="CON_", proteinMatrix$prot)  
                    & !grepl(pattern =".1/sp", proteinMatrix$prot) & !grepl(pattern =";", proteinMatrix$prot), !grepl(pattern ="QC", names(proteinMatrix))]
pM <- pM[grepl("1/sp",pM$prot) | grepl("1/tr", pM$prot),]

#rename proteinMatrix
library(stringr)
nprot <- sapply(pM$prot, function(x){str_split(x, "\\|")[[1]][2]})
ngene <- sapply(pM$prot, function(x){str_split(x, "\\|")[[1]][3]})

pM2 <- data.frame(pM[,2:ncol(pM)],row.names = paste0(nprot,"_",ngene))

t_n <- sapply(names(pM)[2:ncol(pM)], function(x){str_split(x, "_")[[1]][5]})
t_n <- as.character(sapply((t_n), function(x){str_split(x, "\\.")[[1]][1]}))
t_rep <- sapply(names(pM)[2:ncol(pM)], function(x){str_split(x, "_")[[1]][6]})
t_rep <- sapply((t_rep), function(x){str_split(x, "\\.")[[1]][1]})

names(pM2) <- paste0(t_n,"_",t_rep)

#find rep
rep_index <- grep('rep', t_rep)

pM2_rep <- data.frame(NA)
for (i in rep_index){
  pM2_rep <- cbind(pM2_rep,pM2[,i])
}

rep_name <- c("NA")
for (i in rep_index){
  rep_name <- cbind(rep_name,paste0(t_n[i],"_",t_rep[i]))
}
names(pM2_rep) <- rep_name
pM2_rep <- pM2_rep[,-1]

t_rep_name <- c("NA")
for (i in 2:23){
  a <- str_split(rep_name[i], "_")[[1]][1]
  t_rep_name <- cbind(t_rep_name, paste0(a,"_","NA"))
}
for (i in t_rep_name){
  pM2_rep <- cbind(pM2_rep,pM2[,which(names(pM2)==i)])
}

names(pM2_rep)[23:44] <- t_rep_name[2:23]


#### rep取平均值
rep_name <- rep_name[2:23]
for (i in rep_name) {
  pM2 <- pM2[,-which(names(pM2)==i)]
}

average <- data.frame()
for (i in 1:nrow(pM2_rep)){
  for (j in 1:22){
    if( !is.na(pM2_rep[i,j]) & !is.na(pM2_rep[i,j+22])){
      average[i,j] = (pM2_rep[i,j] + pM2_rep[i,j+22]) /2
    }else if(is.na(pM2_rep[i,j]) & !is.na(pM2_rep[i,j+22])){
      average[i,j] = pM2_rep[i,j+22]
    }else if(!is.na(pM2_rep[i,j]) & is.na(pM2_rep[i,j+22])){
      average[i,j] = pM2_rep[i,j]
    }else{
      average[i,j] = NA
    }
  }
}
names(average) <- t_rep_name[2:23]
row.names(average) = paste0(ngene)
t_rep_name <- t_rep_name[2:23]

for(i in t_rep_name){
  pM2[,which(names(pM2)==i)] <- average[,which(names(average)==i)]
}
c <- sapply(names(pM2)[1:ncol(pM2)], function(x){str_split(x, "_")[[1]][1]})
names(pM2) <- c

# type
library(readxl)
info <- read_excel("sample_information_20200507.xlsx")
info2 <- read_excel("patient_informmation_CVDSBB_transfer_20200520_SAA_CRP_3_48h(2)(1).xlsx")
info$label <- info2$Label[match(info$SampleID,info2$SampleID)]
info$clinical <- info2$Clinical_type[match(info$SampleID,info2$SampleID)]


pM3 <- pM2
#NA 填充最小值
pM3[is.na(pM3)] <- min(pM3[,1:348], na.rm = T)
# write.csv(pM3, "pm.csv")

df <- read_csv("skylinenorep.csv")
df <- df[which(df$Protein == "ALBU" | df$Protein == "CRP"),]
df <- df[-2,]
df[4:350] <- log2(df[4:350])
df.albu <- df[2,4:350]
df.crp <- df[1,4:350]

pM4 <- pM3
pM4 <- pM4[, -which(names(pM4) %in% c("52","113"))]
#pM4["P02768_ALBU_HUMAN",] <- df.albu[match(names(pM4),names(df.albu))]
pM4["P02741_CRP_HUMAN",] <- df.crp[match(names(pM4),names(df.crp))]

label_timetag <- info$label[match(names(pM4), info$sample_N)]
clinicallabel <- info$clinical[match(names(pM4), info$sample_N)]

pM4 <- t(pM4)
pM4 <- as.data.frame(pM4)

pM4$label <- label_timetag
pM4$label <- as.factor(pM4$label)
pM4$clinical <- clinicallabel
pM4$clinical <- as.factor(pM4$clinical)
#write.csv(pM4, "pm.csv")
info3 <- read_xlsx("all_prots_p_20200520.xlsx",col_names = F)
source("D:/datamining_library_ge20200306.R")
nm.pm4.2 <- ge.split(names(pM4),"_",2)
nm.pm4.3 <- ge.split(names(pM4),"_",3)
nm.pm4 <- paste0(nm.pm4.2,"_",nm.pm4.3)
nm.pm4[338:339] <- c("label", "clinical")
pM5 <- pM4[, which(nm.pm4 %in% c(info3$...1, "label", "clinical"))]
pM5 <- as.data.frame(pM5)

nm.pm5 <- ge.split(names(pM5),"_",1)
genenm <- read_xlsx("uniprot_gene2.xlsx")
names(pM5) <- genenm$names[match(nm.pm5,genenm$yourlist)]
names(pM5)[49:50] <- c("label", "clinical")


#meandata
# meandata <- data.frame(x= (1:69), y=(1:69))
# 
# for (i in 1:5) {
#   pM4.severe.sub <- pM4.severe[which(as.character(pM4.severe$label) == i),1:69]
#   pM4.severe.sub <- as.data.frame(pM4.severe.sub)
#   a <- c()
#   for (j in 1:69){
#     a[j] <- mean(as.numeric(pM4.severe.sub[,j]))
#   }
#   meandata <- cbind(meandata,a)
#   names(meandata)[i+2] <- paste0("Severe_",i)
# }
# 
# pM4.nonsevere <- pM4[which(pM4$clinical == "Non-Severe"),]
# pM4.nonsevere <- apply(pM4.nonsevere,
#                        2,
#                        as.numeric) %>%
#   data.frame(stringsAsFactors = F)
# for (i in 1:5) {
#   pM4.nonsevere.sub <- pM4.nonsevere[which(pM4.nonsevere$label == i),1:69]
#   a <- apply(pM4.nonsevere.sub,2, mean)
#   meandata <- cbind(meandata,a)
#   names(meandata)[i+7] <- paste0("non-Severe_",i)
# }
# 
# meandata <- meandata[,3:12]
# row.names(meandata) <- names(pM4)[1:69]
#mediandata
pM4.severe <-  pM5[which(pM5$clinical == "Severe"),]
# aov_p = aov(pM4.severe$P02741_CRP~pM4.severe$label)
# p.severe = summary(aov_p)[[1]]["Pr(>F)"][[1]][1]
mediandata <- data.frame(x= (1:48), y=(1:48))

for (i in 1:5) {
  pM4.severe.sub <- pM4.severe[which(as.character(pM4.severe$label) == i),1:48]
  pM4.severe.sub <- as.data.frame(pM4.severe.sub)
  a <- c()
  for (j in 1:48){
    a[j] <- median(as.numeric(pM4.severe.sub[,j]))
  }
  mediandata <- cbind(mediandata,a)
  names(mediandata)[i+2] <- paste0("Severe_",i)
}

pM4.nonsevere <- pM5[which(pM5$clinical == "Non-Severe"),]
# aov_p = aov(pM4.nonsevere$P02741_CRP~pM4.nonsevere$label)
# p.nonsevere = summary(aov_p)[[1]]["Pr(>F)"][[1]][1]
# pM4.nonsevere <- apply(pM4.nonsevere,
#                        2,as.numeric) %>%data.frame(stringsAsFactors = F)
for (i in 1:5) {
  pM4.nonsevere.sub <- pM4.nonsevere[which(pM4.nonsevere$label == i),1:48]
  a <- apply(pM4.nonsevere.sub,2, median)
  mediandata <- cbind(mediandata,a)
  names(mediandata)[i+7] <- paste0("non-Severe_",i)
}

pM4.nonCOVID <- pM5[which(pM5$clinical == "non-COVID"),]
mediandata <- cbind(mediandata,apply(pM4.nonCOVID[,1:48],2, median))
names(mediandata)[13] <- paste0("non-COVID")
pM4.Healthy <- pM5[which(pM5$clinical == "Healthy"),]
mediandata <- cbind(mediandata,apply(pM4.Healthy[,1:48],2, median))
names(mediandata)[14] <- paste0("Healthy")


mediandata <- mediandata[,c(14,13,8:12,3:7)]

# mediandata <- mediandata[-which(row.names(mediandata) %in%c("IGF2_HUMAN","NOX5_HUMAN", "SIA8A_HUMAN","FHR2_HUMAN")), ] 
# # median.severe <- mediandata[,grep("^Severe", names(mediandata))]
# mediandata_del <- mediandata[-which(row.names(mediandata) %in% c("NOX5_HUMAN", "ACD_HUMAN","K0513_HUMAN", "YAED1_HUMAN", "ANT3_HUMAN", 
#                                                                 "CO8G_HUMAN","CO8A_HUMAN","LCAT_HUMAN","C1R_HUMAN","FA78B_HUMAN",
#                                                                 "MRC1_HUMAN","RM45_HUMAN","CO2A1_HUMAN","ZN672_HUMAN")), ] 
# meandata_del <- meandata[-which(row.names(meandata)%in% c("NOX5_HUMAN", "ACD_HUMAN","K0513_HUMAN", "YAED1_HUMAN", "ANT3_HUMAN", 
#                                                           "CO8G_HUMAN","CO8A_HUMAN","LCAT_HUMAN","C1R_HUMAN","FA78B_HUMAN",
#                                                       


# rename
# infoname <- read.csv("sp_protein_gene.csv")
# infoname$id2 <- sapply(infoname$id, function(x){str_split(x, "\\|")[[1]][2]})
# infoname$id3 <- sapply(infoname$id, function(x){str_split(x, "\\|")[[1]][3]})
# median_name <- infoname$id2[match(row.names(mediandata_del), infoname$id3)]
# median_gene <- as.character(infoname$gene[match(row.names(mediandata_del), infoname$id3)])  
# row.names(mediandata_del) <- paste0(median_name,"_",median_gene)
# mediandata_del <- mediandata_del[-which(row.names(mediandata_del)%in% c("P02768_ALB","P02741_CRP")),]
# #draw
# write.csv(mediandata_del, "heatmap_originalname.csv")
# mediandata_del.score <- apply(mediandata_del,2, scale)
# write.csv(mediandata_del.score, "heatmap_scale.csv")

# pheatmap(mediandata, color = c(brewer.pal(9,"YlOrRd")[1:6]), fontsize_col = 8,
#          # annotation_col = pickann_col,
#          # annotation_colors = pickann_colors, 
#          scale = "row",
#          cluster_rows = T, cluster_cols = F,show_rownames = T, show_colnames = T, 
#          filename = "heatmap_median.pdf",main = "heatmap",width=8,height=12)

pheatmap(mediandata, color = c(brewer.pal(11,"RdYlBu")[9:7],"azure1",brewer.pal(11,"RdYlBu")[4:1]), fontsize_col = 12,
         # annotation_col = pickann_col,
         # annotation_colors = pickann_colors, 
         scale = "row",
         cluster_rows = T, cluster_cols = F,show_rownames = T, show_colnames = T, 
         filename = "heatmap_medianv5.pdf",main = "heatmap",width=8,height=10)
