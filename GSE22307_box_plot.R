
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))
install.packages(c("dplyr", "tidyverse", "ggplot2", "ggpubr", "stringr"))

library(GEOquery)
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)

# data retival 

# matrix retival
gse<-getGEO("GSE22307", GSEMatrix = T)
eSet <- gse[[1]]
raw_data<-exprs(eSet)

normalize_data<-normalizeBetweenArrays(raw_data, method = "quantile")
normalize_data1<-log2(normalize_data+1)

# probe annotation data
feature_data<-fData(eSet)

# change probes name to symbols
probe_ID<-feature_data$ID
geneSymbol<-feature_data$`Gene Symbol`

# Replace rownames of normalized_data with corresponding gene symbols
rownames(normalize_data1) <- geneSymbol[match(rownames(normalize_data1), probe_ID)]

# phenotype data
pheno_data<-pData(eSet)
col_data<-pheno_data[,c("title","time point (day):ch1","geo_accession")]
col_data<-col_data%>%mutate(group = ifelse(`time point (day):ch1`==0, "0 Days of DSS treatment", paste0( `time point (day):ch1`, " Days of DSS treatment")))
col_data$day_of_DSS_treatment<-as.numeric(col_data$`time point (day):ch1`)

# OC for quntile normalization process
color_pal<- c("green", "pink","red3","red4")
boxplot(log2(raw_data+1), col= color_pal[as.factor(col_data$day_of_DSS_treatment)]) 
boxplot(log2(normalize_data+1), col= color_pal[as.factor(col_data$day_of_DSS_treatment)])

# calculate Z socre
norm_z_score<- t(scale(t(normalize_data1),center =T, scale =T ))
norm_z_score<-as.data.frame(na.omit(pmax(pmin(norm_z_score, 2), -2)))
norm_z_score$gene<-rownames(norm_z_score)

# long data
data_long<-norm_z_score %>%
  gather(key ='samples', value = 'Z_score',-gene)%>%
  left_join(.,col_data, by=c("samples"="geo_accession"))

#box plots function 

boxplot_gene_expression <- function(gene_name) {
  p <- data_long %>%
    filter(gene == gene_name) %>%
    ggplot(aes(x = group, y = Z_score, fill = group)) +
    geom_boxplot() +
    geom_jitter(shape = 20, position = position_jitter(0.1)) +  # Adding jitter
    stat_compare_means(
      comparisons = list(
        c("0 Days of DSS treatment", "2 Days of DSS treatment"),
        c("0 Days of DSS treatment", "4 Days of DSS treatment"),
        c("0 Days of DSS treatment", "6 Days of DSS treatment")
      ),
      method = "t.test",
      label = "p",
      bracket.size = 0.9,
      vjust = -0.5,
      hide.ns = TRUE
    ) +
    theme_classic() +
    scale_fill_manual(values = c(
      "0 Days of DSS treatment" = "green",
      "2 Days of DSS treatment" = "pink",
      "4 Days of DSS treatment" = "red3",
      "6 Days of DSS treatment" = "red4"
    )) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 90,size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold")
    ) +
    labs(title = gene_name)
  
  return(p)
}

boxplot_gene_expression("S100a9")
boxplot_gene_expression("S100a8")
boxplot_gene_expression("Sox9")
