#!/usr/bin/Rscript
## Author: Merai Dandouch
## meraid@bu.edu
## BU BF528
## Project 3

library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('ggplot2')
library('ggVennDiagram')

#### load count data ####
#' @param filename full file path as a string of the counts file
#' @param sample_name sample for the SRR file, ex: SRR1177966
#' 
#' @return A _data frame_ with gene ids as row names. A tibble will **not** work 
#' @details As always, we need to load our data and start to shape it into the 
#' form we need for our analysis. Selects only the columns named "gene", "vP0_1", 
#' "vP0_2", "vAd_1", and "vAd_2" from the counts file. 
load_countdata <- function(filename, sample_name) {
  count_data <- read.table(filename, sep = "\t", skip = 1, header = T) #read.table the samples 
  names(count_data)[7] <- sample_name #change count column name to sample name 
  count_data %>% dplyr::select(Geneid, sample_name) -> count_data #select Geneid and count columns only 
  return (count_data)
} 

#### merge count data ####
#' @param count_list list of count files from
#' 
#' @return _counts_,  a merged _data frame_ with gene ids as row names and samples 
#' as column names. A tibble will **not** work 
merge_countdata <- function(count_list){
  counts <- Reduce(function(x, y) merge(x, y, all=TRUE), count_list) 
  return (counts)
}

#### remove outliers from boxplot results ####
#' @param count_matrix matrix of gene ids and sample counts 
#' 
#' @return _results_,  a _data frame_ with gene ids as row names and samples 
#' as column names where outliers have been removed 
#' @details The boxplot() object has stat information where one can check and see
#' which values are the outliers. I removed these values from count_matrix. However,
#' this caused some rows to contain NANs so later, I removed the rows with NANs.
remove_outliers <- function(count_matrix) {
  data <- pivot_longer(count_matrix, cols = starts_with('SRR'), names_to = 'Samples', values_to = 'Counts')
  out <- boxplot.stats(data$Counts)$out
  out_ind <- which(data$Counts %in% c(out))
  results <- data[-out_ind, ]
  results <- pivot_wider(results, names_from = 'Samples', values_from = 'Counts')
  return(results)
}

#### box plot figure ####
#' @param count_matrix matrix of gene ids and sample counts where outliers have 
#' been removed
#' 
#' @return _p_,  a gg plot object that contain a boxplot figure 
make_boxplots <- function(count_matrix){
  data <- pivot_longer(count_matrix, cols = starts_with('SRR'), names_to = 'Samples', values_to = 'Counts')
  p <- ggplot(data, aes(x=Samples, y=Counts, color=Samples)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90))
  return(p)
}

#' Filter out genes with zero variance
#' @param verse_counts matrix of gene ids and sample counts where outliers have 
#' been removed
#'
#' @return tibble: a (n x m) tibble of raw reads with genes that have 
#' zero variance across samples removed
#' 
#' @note (g >= n)
#' 
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`
filter_zero_var_genes <- function(verse_counts) {
  nonzero_genes <- apply(verse_counts[-1], 1, var) #apply var on every row but exclude the first col
  filtered_counts <- verse_counts[nonzero_genes!=0,]
  filtered_counts <- subset(filtered_counts, rowSums(filtered_counts==0)==0)
  filtered_counts
  return(filtered_counts)
}

#' Perform a DESeq2 analysis of rna seq data
#'
#' @param counts matrix of gene ids and sample counts where outliers have 
#' been removed
#' @param info The coldata variable describing the experiment, a dataframe.
#' @param exp The vehicle introduced by chemical causing environmental stress 
#' response in organism.
#'
#' @return A dataframe of DESeq results. It has a header describing the 
#' stats (padj, log2FoldChange, pvalue...) of results with genes as row names. 
#' @details This function is based on the DESeq2 User's Guide. These links describe 
#' the inputs and process we are working with. The output we are looking for comes 
#' from the DESeq2::results() function.
#' https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
#' https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
#'
#' @examples run_deseq(master_matrix, info, "AhR")
run_deseq <- function(counts, info, exp){
  counts <- as_tibble(counts) %>% column_to_rownames(var = "Geneid") %>% as.matrix()
  # create the DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = info,
    design= ~ mode_of_action)
  
  # relevel mode_of_action as factor
  dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')
  
  # run DESeq
  dds <- DESeq(dds)
  res <- results(dds, contrast=c('mode_of_action', paste(info[2,2]),'Control'))
  res <- lfcShrink(dds, coef=2)
  
  #Save Normalized Counts
  write.csv(counts(dds,normalized=TRUE), paste("data/deseq_norm_", exp,".csv",sep=""))
  return(res)
}

#### histogram figure ####
#' @param results A dataframe of DESeq results. It has a header describing the 
#' stats (padj, log2FoldChange, pvalue...) of results with genes as row names. 
#' 
#' @return _p_,  a gg plot object that contain a histogram figure 
#' @details Please Note: the results dataframe being passed through the parameter
#' has already been filtered for padj < 0.05 and arranged by ascending pvalues
create_histogram <- function(results) {
  results %>% filter(padj < 0.05) 
  p <- ggplot(results, mapping = aes(x = log2FoldChange)) +
    geom_histogram(colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") + 
    geom_vline(aes(xintercept=mean(log2FoldChange)),
               color="blue", linetype="dashed", size=1) + 
    labs(title = "Tox Group 2 Fold Change Histogram", x = "Log2 Fold Change", y = "Frequency")
  print(p)
  return(p)
}

#### scatterplot figure ####
#' @param results A dataframe of DESeq results. It has a header describing the 
#' stats (padj, log2FoldChange, pvalue...) of results with genes as row names. 
#' 
#' @return _p_,  a gg plot object that contain a scatterplot figure 
#' @details Please Note: the results dataframe being passed through the parameter
#' has already been filtered for padj < 0.05 and arranged by ascending pvalues.
#' As you can see all the points the figure are above the significance level, or
#' in this case the horizontal line.
create_scatterplot <- function(results) {
  p <- ggplot(results, aes(x = log2FoldChange, y = -log(pvalue)), color = pvalue < 0.05) + 
    geom_point() + 
    geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
    scale_color_manual(values=c("black", "red")) + 
    ggtitle("Fold Change Scatter Plot") + 
    theme_classic()
  return(p)
}


#Read in count files prodcued by featureCounts() from subread package 
count_list <- vector(mode = "list", length = 9) # 9 samples 
count_list[[1]] <- load_countdata('feature_counts/SRR1177966Aligned.sortedByCoord.out.bam.txt', 'SRR1177966')
count_list[[2]] <- load_countdata('feature_counts/SRR1177969Aligned.sortedByCoord.out.bam.txt', 'SRR1177969')
count_list[[3]] <- load_countdata('feature_counts/SRR1177970Aligned.sortedByCoord.out.bam.txt', 'SRR1177970')
count_list[[4]] <- load_countdata('feature_counts/SRR1177993Aligned.sortedByCoord.out.bam.txt', 'SRR1177993')
count_list[[5]] <- load_countdata('feature_counts/SRR1177994Aligned.sortedByCoord.out.bam.txt', 'SRR1177994')
count_list[[6]] <- load_countdata('feature_counts/SRR1177995Aligned.sortedByCoord.out.bam.txt', 'SRR1177995')
count_list[[7]] <- load_countdata('feature_counts/SRR1177998Aligned.sortedByCoord.out.bam.txt', 'SRR1177998')
count_list[[8]] <- load_countdata('feature_counts/SRR1178001Aligned.sortedByCoord.out.bam.txt', 'SRR1178001')
count_list[[9]] <- load_countdata('feature_counts/SRR1178003Aligned.sortedByCoord.out.bam.txt', 'SRR1178003')

# Treatment Manipulation: Merge all treatment samples list into one dataframe, extract treatment names, and filter
treatment_matrix <- merge_countdata(count_list) %>% tibble::as_tibble()
treatment_samples <- colnames(treatment_matrix[-1]) #as_tibble() producing extra col
treatment_matrix <- filter_zero_var_genes(treatment_matrix) %>% remove_outliers() 
treatment_matrix <- na.omit(treatment_matrix)
bp <- make_boxplots(treatment_matrix)
ggsave('data/bp.png', plot = bp)

# Control Manipulation: Read metadata, subset control samples from metadata and pull only needed for tox group 2 from control_counts.csv
metadata <- read.csv("/project/bf528/project_3/groups/group_2_rna_info.csv")
write.csv(metadata, "data/metadata.csv")
control_samples <- metadata[metadata$mode_of_action == 'Control', ] %>% dplyr::pull(Run)
control_matrix <- read.csv("/project/bf528/project_3/samples/control_counts.csv") %>% dplyr::select(Geneid, control_samples)

# Merge Treatment and Control Matrix by cols
master_matrix <- merge(treatment_matrix, control_matrix, by="Geneid") # first 9 cols are treatment cols merged with next 9 control cols 

# Subset master matrix to include only grouped samples for a given chemical/vehicle experiment before running DESeq2. 
ahr <- master_matrix %>% dplyr::select(Geneid, as.character(metadata[metadata$vehicle == 'CMC_.5_%',] %>% dplyr::pull(Run)))
ahr_info <- metadata[metadata$vehicle == 'CMC_.5_%',] #pull metadata for AhR vehicle
rownames(ahr_info) <- NULL 
ahr_res <- run_deseq(ahr, ahr_info, 'AhR')
write.csv(ahr_res, "data/ahr_res.csv")

car <- master_matrix %>% dplyr::select(Geneid, as.character(metadata[metadata$vehicle == 'CORN_OIL_100_%',] %>% dplyr::pull(Run)))
car_info <- metadata[metadata$vehicle == 'CORN_OIL_100_%',] #pull metadata for CAR vehicle
rownames(car_info) <- NULL 
car_res <- run_deseq(car, car_info, 'CAR')
write.csv(car_res, "data/car_res.csv")

cyto <- master_matrix %>% dplyr::select(Geneid, as.character(metadata[metadata$vehicle == 'SALINE_100_%',] %>% dplyr::pull(Run)))
cyto_info <- metadata[metadata$vehicle == 'SALINE_100_%',] #pull metadata for cytotoxic vehicle
rownames(cyto_info) <- NULL 
cyto_res <- run_deseq(cyto, cyto_info, 'Cyto')
write.csv(cyto_res, "data/cyto_res.csv")

#Find signficant genes and create histogram, scatterplot for 3 groups 
ahr_res %>% as.data.frame() %>% arrange(pvalue) %>% dplyr::filter(padj < 0.05) -> ahr_res
ahr_res %>% dplyr::slice(1:10) %>% cbind(Group = rep("BETA-NAPHTHOFLAVONE",length(10)))-> AHR_top10DE
p1 <- create_histogram(ahr_res) 
p2 <- create_scatterplot(ahr_res)
ggsave("data/ahr_histogram.png", plot = p1)
ggsave("data/ahr_scatterplot.png", plot = p2)

car_res %>% as.data.frame() %>% arrange(pvalue) %>% dplyr::filter(padj < 0.05)-> car_res
car_res %>% dplyr::slice(1:10) %>% cbind(Group = rep("ECONAZOLE",length(10)))-> CAR_top10DE
p3 <- create_histogram(car_res) 
p4 <- create_scatterplot(car_res)
ggsave("data/car_histogram.png", plot = p3)
ggsave("data/car_scatterplot.png", plot = p4)

cyto_res %>% as.data.frame() %>% arrange(pvalue) %>% dplyr::filter(padj < 0.05) -> cyto_res
cyto_res %>% dplyr::slice(1:10) %>% cbind(Group = rep("THIOACETAMIDE",length(10)))-> CYTO_top10DE
p5 <- create_histogram(cyto_res) 
p6 <- create_scatterplot(cyto_res)
ggsave("data/cyto_histogram.png", plot = p5)
ggsave("data/cyto_scatterplot.png", plot = p6)

#save top 10 differentially expressed genes for the 3 groups (AhR, CAR, Cytotoxic)
rbind(AHR_top10DE, CAR_top10DE, CYTO_top10DE) %>% write.csv("data/top10DE.csv")

