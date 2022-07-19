library(tidyr)
library(tidyverse)
library(dplyr)
library(limma)
library(ggplot2)
library(glue)
library(gridExtra)
library(ggpubr)

#-------------------------------------Section 5--------------------------------------------------#
# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_2_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

run_limma <- function(chemical, samples, rma){
  # subset the full expression matrix to just those in this comparison
  rma.subset <- rma[paste0('X',samples$array_id[samples$chemical == 'Control'|samples$chemical==chemical])]

  # construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(
    ~factor(
      samples$chemical,
      levels=c('Control',chemical)
    )
  )
  colnames(design) <- c('Intercept', chemical)
  
  # run limma
  fit <- lmFit(rma.subset, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH', sort.by='P')
  #filter t
  t_filter <- filter(t, adj.P.Val<0.05)
  top10 <- t_filter[1:10,]
  
  #write out the results to file
  write.csv(t_filter,glue('{chemical}_limma_results.csv'))
  write.csv(top10,glue('top10_{chemical}_limma_results.csv'))

  plot_limma(t,t_filter,top10, chemical)
}


plot_limma <- function(t, t_filter, top10, chemical){
  #histogram
  plot_h <- t_filter %>%
    ggplot(aes(logFC)) +
    geom_histogram( binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    geom_vline(aes(xintercept = 0),
               color="blue", linetype="dashed", size=1) +
    geom_density(alpha=.2, fill="#FF6666") +
    labs(x=glue("Log2 Fold Change (Control vs. {chemical})"),
         y="Frequency") 
  
  plot_s <- t%>%
    ggplot(aes(x=logFC, y=-log10(P.Value), color=P.Value<0.05)) +
    geom_point() +
    scale_color_manual(values=c("black", "red")) +
    labs(x=glue("Log2 Fold Change (Control vs. {chemical})"),
         y="Normal p-value")
  return(grid.arrange(grobs = list(plot_h,plot_s)))
}

for (chemical in c('BETA-NAPHTHOFLAVONE', 'ECONAZOLE', 'THIOACETAMIDE')){
  run_limma(chemical, samples, rma)
}
run_limma('THIOACETAMIDE', samples, rma)

#----------------------------------------------------------------------------------------



#-------------------------------------Section 6--------------------------------------------------#
###########read all required files
map <- read.csv('/project/bf528/project_3/refseq_affy_map.csv')

#read deseq files and merge it with map file by 'refseq'
read_deseq <- function(path, map){
  result <- read.csv(glue("/projectnb/bf528/users/hedgehog_2022/project3/Programmer/project-3-hedgehog-1/data/{path}_res.csv")) %>%
    rename(REFSEQ = X) %>%
    merge(select(map, c(REFSEQ, SYMBOL)), by='REFSEQ')
  return(result)
}
beta_deseq <- read_deseq('ahr', map)
econ_deseq <- read_deseq('car', map)
thioa_deseq <- read_deseq('cyto', map)

#read limma files and merge it with map file by 'probeid'
read_limma <- function(path, map){
  result <- read.csv(glue("/projectnb/bf528/users/hedgehog_2022/project3/Analyst/{path}_limma_results.csv")) %>%
    rename(PROBEID = X)%>%
    merge(select(map, c(PROBEID, SYMBOL)), by='PROBEID')
  return(result)
}
beta_limma <- read_limma('BETA-NAPHTHOFLAVONE', map)
econ_limma <- read_limma('ECONAZOLE', map)
thioa_limma <- read_limma('THIOACETAMIDE', map)



################6.3####################
#calculate concordance
concordance <- function(file_deseq,file_limma){
  n0 <- length(intersect(file_deseq$SYMBOL, file_limma$SYMBOL))
  n1 <- nrow(file_deseq)
  n2 <- nrow(file_limma)
  #N=n1+n2-n0
  # x=0
  # intersect <- x+(n1-x)*(n2-x)/(N-x)
  # concordance <- 2*intersect/(n1+n2)
  N=54879
  x <- (N*n0-n1*n2)/(n0+N-n1-n2)
  concordance <- (2*x)/(n1+n2)
  
  return(concordance)
}

#get result for three analyses
conc_beta <- concordance(beta_deseq, beta_limma) 
conc_econ <- concordance(econ_deseq, econ_limma) 
conc_thioa <- concordance(thioa_deseq, thioa_limma) 




###################6.4#####################
#plot concordance vs number of DE genes for both analysis
plot_overall <- function(table, type){
  ggplot(table, aes(Treatment, Concordance)) +
    geom_point() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_text(label=c('BETA-NAPHTHOFLAVONE', 'ECONAZOLE', 'THIOACETAMIDE')) +
    geom_smooth(method = "lm") +
    stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
    labs(x=glue('Treat Effect (number of DEGs from {type})'), y="Concordance of DEG") %>%
    return()
}

#deseq
deseq_table <- data.frame("Concordance"=c(conc_beta, conc_econ, conc_thioa), 
                          "Treatment"=c(nrow(beta_deseq), nrow(econ_deseq), nrow(thioa_deseq)),
                          "Chemical"=c('BETA-NAPHTHOFLAVONE', 'ECONAZOLE', 'THIOACETAMIDE'))
plot_deseq <- plot_overall(deseq_table, "RNA-seq")
#limma
limma_table <- data.frame("Concordance"=c(conc_beta, conc_econ, conc_thioa), 
                          "Treatment"=c(nrow(beta_limma), nrow(econ_limma), nrow(thioa_limma)),
                          "Chemical"=c('BETA-NAPHTHOFLAVONE', 'ECONAZOLE', 'THIOACETAMIDE'))
plot_limma <- plot_overall(limma_table, "microarray")




##################6.5#####################
######split beta deseq file
#deseq
beta_deseq_mean<-median(beta_deseq$baseMean)
above_beta_deseq <- subset(beta_deseq, baseMean >= beta_deseq_mean) 
below_beta_deseq <- subset(beta_deseq, baseMean < beta_deseq_mean)
#limma
beta_limma_mean<-median(beta_limma$AveExpr)
above_beta_limma <- subset(beta_limma, AveExpr >= beta_limma_mean) 
below_beta_limma <- subset(beta_limma, AveExpr < beta_limma_mean)

######split econ deseq file
#seq
econ_deseq_mean<-median(econ_deseq$baseMean)
above_econ_deseq <- subset(econ_deseq, baseMean >= econ_deseq_mean) 
below_econ_deseq <- subset(econ_deseq, baseMean < econ_deseq_mean)
#limma
econ_limma_mean<-median(econ_limma$AveExpr)
above_econ_limma <- subset(econ_limma, AveExpr >= econ_limma_mean) 
below_econ_limma <- subset(econ_limma, AveExpr < econ_limma_mean)

######split thioa deseq file
#deseq
thioa_deseq_mean<-median(thioa_deseq$baseMean)
above_thioa_deseq <- subset(thioa_deseq, baseMean >= thioa_deseq_mean) 
below_thioa_deseq <- subset(thioa_deseq, baseMean < thioa_deseq_mean)
#limma
thioa_limma_mean<-median(thioa_limma$AveExpr)
above_thioa_limma <- subset(thioa_limma, AveExpr >= thioa_limma_mean) 
below_thioa_limma <- subset(thioa_limma, AveExpr < thioa_limma_mean)

##########calculate sub concordance
#beta concoredance
conc_beta_above <- concordance(above_beta_deseq, above_beta_limma) #0.0139561
conc_beta_below <- concordance(below_beta_deseq, below_beta_limma) #0.01397574
#econ concordance
conc_econ_above <- concordance(above_econ_deseq, above_econ_limma) #0.1936766
conc_econ_below <- concordance(below_econ_deseq, below_econ_limma) #0.1927679

#thioa concordance
conc_thioa_above <- concordance(above_thioa_deseq, above_thioa_limma) #0.4279785
conc_thioa_below <- concordance(below_thioa_deseq, below_thioa_limma) #0.4204899




################6.6#######################
#concordance for three analysis
#          beta               econ                thioa
# overall   conc_beta          conc_econ           conc_thioa
# above     conc_beta_above    conc_econ_above     conc_thioa_above
# below     conc_beta_below    conc_econ_below     conc_thioa_below

bar_plot <- data.frame("Concordance"=c(conc_beta,conc_beta_above,conc_beta_below,conc_econ,conc_econ_above,conc_econ_below,conc_thioa,conc_thioa_above,conc_thioa_below),
                       "Level"=rep(c('Overall', 'Above', 'Below'), 3),
                       "Chemical"= c(rep("BETA-NAPHTHOFLAVONE" , 3) , rep("ECONAZOLE" , 3) , rep("THIOACETAMIDE" , 3))) %>%
  ggplot(aes(fill=Level, y=Concordance, x=Chemical)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_text(aes(label = round(Concordance, 3)), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.position = c(.10, .80)) +
  labs(x='Treatment')




