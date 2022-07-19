library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tidyverse)
library(gplots)
library(stats)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(grid)

#Reading the deseq norm counts

ahr<-read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Programmer/project-3-hedgehog-1/data/deseq_norm_AhR.csv")
car<-read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Programmer/project-3-hedgehog-1/data/deseq_norm_CAR.csv")
cyto<-read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Programmer/project-3-hedgehog-1/data/deseq_norm_Cyto.csv")

#merging the columns
norm<-cbind(ahr,car,cyto)
row.names(norm)<-norm[[1]]
#tidying the data
norm<-norm[,-c(1,8,15)]
#filtering
#mean>20
Mean<-apply(norm,1,mean)
Filter1 <- filter(norm, Mean>20)

#Chi squared test
df = length(Filter1)-1 #degree of freedom = n-1
chilower = qchisq((0.01)/2, df) #pvalue = 0.01, lower tail
chiupper = qchisq((1 - 0.01)/2, df, lower.tail = FALSE) #upper tail
Variance<-apply(Filter1, 1, var) #Variance of genes
Test_statistic <- (df*Variance/median(Variance))
Filter2 <-filter(Filter1, Test_statistic > chiupper | Test_statistic < chilower) 

#Coefficient of Variation
CV<-function(x){
  sd(x)/mean(x)
}
coeff_var<-apply(Filter1,1,CV)
Filter3<-dplyr::filter(Filter1, coeff_var > 0.2)

#Heatmap
Plot<-data.matrix(Filter3)
colors = brewer.pal(n = 11, name = "PiYG")
colors = colorRampPalette(colors)(50)
colors = rev(colors)
metadata <- data.frame(Sample=c("AhR-BETA-NAPHTHOFLAVONE", "AhR-BETA-NAPHTHOFLAVONE", "AhR-BETA-NAPHTHOFLAVONE", "Control", "Control", "Control",
                                "CAR/PXR-ECONAZOLE", "CAR/PXR-ECONAZOLE", "CAR/PXR-ECONAZOLE", "Control", "Control", "Control",
                                "Cytotoxic-THIOACETAMIDE", "Cytotoxic-THIOACETAMIDE", "Cytotoxic-THIOACETAMIDE","Control", "Control", "Control")
)
rownames(metadata) <- colnames(Plot)





## Create the heatmap:
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=1, name="vp", just=c("right","top"))), action="prepend")
heatmap <- pheatmap(Plot, scale = "row", annotation_col = metadata, color = colors,fontsize_row = 4,border_color = NA,
                    clustering_distance_cols="euclidean", main = "Clustering of Treatment and Control Samples based on MOA ", show_rownames = FALSE)

setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=10))
grid.text("Genes", x=-0.07, rot=90, gp=gpar(fontsize=10))

#DAVID analysis
ahr_deseq<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Biologist/data/ahr_res.csv")
car_deseq<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Biologist/data/car_res.csv")
cyto_deseq<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Biologist/data/cyto_res.csv")
ahr_limma<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Analyst/BETA-NAPHTHOFLAVONE_limma_results.csv")
car_limma<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Analyst/ECONAZOLE_limma_results.csv")
cyto_limma<- read_csv("/projectnb/bf528/users/hedgehog_2022/project3/Analyst/THIOACETAMIDE_limma_results.csv")
#function to filter genes -- pvalue<0.05 and logfoldchange> log2(1.5)
filter_func <- function (filename) {
  return (filter(filename, padj<0.05 & abs(log2FoldChange)>log2(1.5)))
}
filter_function <- function (filename) {
  return (filter(filename, adj.P.Val<0.05 & abs(logFC)>log2(1.5)))
}
ahr_deseqfilt<- filter_func(ahr_deseq)
car_deseqfilt<- filter_func(car_deseq)
cyto_deseqfilt<- filter_func(cyto_deseq)
ahr_limmafilt<- filter_function(ahr_limma)
car_limmafilt<- filter_function(car_limma)
cyto_limmafilt<- filter_function(cyto_limma)

# write_csv(ahr_deseqfilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/ahr-d.csv")
# write_csv(ahr_limmafilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/ahr-l.csv")
# write_csv(car_deseqfilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/car-d.csv")
# write_csv(car_limmafilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/car-l.csv")
# write_csv(cyto_deseqfilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/cyto-d.csv")
# write_csv(cyto_limmafilt, "/projectnb/bf528/users/hedgehog_2022/project3/Biologist/cyto-l.csv")




