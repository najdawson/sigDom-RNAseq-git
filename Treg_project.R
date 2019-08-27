############################
#title: "Treg project"
#author: "German Novakovskiy"
#date: "April 27, 2019"
############################

#cleaning the environment
rm(list=ls())

#loading required packages
library(limma) 
library(edgeR)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble)
library(pheatmap)
library(sva)
library(fgsea)
library(stats)
library(gplots)
library(RSvgDevice)

#Loading the data (change folders to apropriate paths)
files <- list.files("/home/german/Nick_Treg_project/Nick_filtered_counts/counts/")
x <- readDGE(files, path = "/home/german/Nick_Treg_project/Nick_filtered_counts/counts/", columns=c(1,2)) 
#removing meta tags
x$counts <- x$counts[1:(nrow(x$counts)-5),]

#Loading metadata and specifying sex factor variable
meta_data <- read.table(file = "~/Nick_Treg_project/SigDom RNAseq data/CAR-Treg RNAseq sample info.csv", sep=",", header=T, row.names = 1)
meta_data$SampleName <- paste(meta_data$Donor, meta_data$Construct, sep=".")
sex <-meta_data$Donor
sex[sex %in% c(9649, 967)] <- "F"
sex[sex %in% c(960)] <- "M"
meta_data$sex <- sex
meta_data$FileName <- sapply(as.character(meta_data$FileName), 
                             function(x){ 
                               paste(unlist(strsplit(x, split="[.]"))[1], "counts", sep=".") 
                               })

x$counts <- x$counts[,meta_data$FileName]
x$samples <- x$samples[meta_data$FileName,]

group <- meta_data$Construct
x$samples$group <- group

#Filtering lowly expressed genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE] #14098    48

#Normalizing for the library size
x <- calcNormFactors(x, method = "TMM")

#For plotting MDS plot (colored for donors)
x$samples$batch <- as.factor(meta_data$Batch)
x$samples$sex <- as.factor(meta_data$sex)
x$samples$stimulus <- as.factor(meta_data$Stimulus)
x$samples$Cell <- as.factor(meta_data$Cell)
x$samples$Donor <- as.factor(meta_data$Donor)

lcpm <- cpm(x, log=TRUE)
col.donor <- x$samples$Donor
levels(col.donor) <-  brewer.pal(nlevels(col.donor), "Set1")
col.donor <- as.character(col.donor)
plotMDS(lcpm, cex=0.75, labels = group, col = col.donor, dim.plot = c(1,2)) #first and second dimension
plotMDS(lcpm, cex=0.75, labels = group, col = col.donor, dim.plot = c(3,4)) #third and fourth

#############################################################################################################
#DGE ANALYSIS
#############################################################################################################
meta_data$Batch <- as.factor(meta_data$Batch)
meta_data$sex <- as.factor(meta_data$sex)
meta_data$Donor <- as.factor(meta_data$Donor)

#Different questions and corresponding groups for analysis
#1 Which transcriptional pathways are differentially activated in CARs that work in vivo vs. those that don’t. Try two separate analyses:
# Analysis A: (stringent analysis)
group1A1 <- c("CD28wt") 
group1A2 <- c("PD1", "TNFR2", "C4wt", "C4mut", "3zeta") 
# Analysis B: (less stringent analysis)

#2 Why do 41BB/TNFR2 destabilize Treg phenotype after CAR stimulation?
#Analysis A: (top performing CARs)
group2A1 <- c("CD28wt") 
group2A2 <- c("41BB", "TNFR2") 

#3 Why does CD28wt stimulate proliferation and suppression but CD28mut does not? 
#Analysis A: (Tregs)
group3A1 <- c("CD28wt") 
group3A2 <- c("CD28mut") 

#4 Why does CD3zeta-only activate CAR stimulation but PD1+CD3zeta does not?
group4A1 <- c("3zeta") 
group4A2 <- c("PD1") 

#5 Which transcriptional pathways are differentially activated in a TCR stim vs a CAR stim:
group5A1 <- c("CD28wt") 
group5A2 <- c("NGFR") 

group5B1 <- c("CD28wt_Tconv") 
group5B2 <- c("NGFR_Tconv") 

#6 Which transcriptional pathways are differentially activated in Treg vs Tconv CAR stimulation? 
group6A1 <- c("CD28wt") 
group6A2 <- c("CD28wt_Tconv") 

#Loading gmt files for enrichment analysis
termGO <- gmtPathways("/home/german/Nick_Treg_project/Human_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt")
termTF <- gmtPathways("/home/german/Nick_Treg_project/RegNetworkStrong.gmt")
termH <- gmtPathways("/home/german/Nick_Treg_project/h.all.v6.2.symbols.gmt")

###############################
#Help functions
###############################

###############
#function for DGE for specified question (groups)
limma_DGE_group_analysis <- function(expr_mat, meta_data, group1, group2){
  #Introducing the corresponding column to meta data 
  md <- meta_data
  md <- md %>% filter(Construct %in% c(group1, group2))
  md$Analysis <- integer(nrow(md))
  md$Analysis[which(md$Construct %in% group1)] <- 1
  md$Analysis[which(md$Construct %in% group2)] <- 2
  md$Analysis <- as.factor(md$Analysis)  
  
  #specifying design matrix; regress out Donor (since it's nested in batch, it will take into account batch effects); 
  #no intercept (0) so easier to look at differences between groups
  designMatrix <- model.matrix(~0 + Analysis + Donor, md)
  #contrast matrix to specify what type of differences we are looking for 
  contrastMatrix <- makeContrasts(
    A1 = Analysis1 - Analysis2,
    levels = colnames(designMatrix))
  
  #subsetting the corresponding samples
  expr_mat <- expr_mat[,md$FileName]
  
  #voom -  estimates the mean-variance relationship and uses this to compute appropriate observation-level weights
  v <- voom(expr_mat, designMatrix, plot=FALSE)
  #Fit linear model for each gene
  vfit <- lmFit(v, designMatrix)
  #fit the contrast using the original fitted model
  vfit <- contrasts.fit(vfit, contrasts=contrastMatrix)
  #apply eBayes() for moderated statistics
  efit <- eBayes(vfit)
  
  #0.05 adjusted p-value
  cutoff <- 5e-02 
  #Multiple Test Correction method by default is BH (equivalent to FDR)
  dt <- decideTests(efit, p.value = cutoff, lfc = 1)
  dt <- as.data.frame(dt)
  #table with all genes and their DGE stats
  DEgenes <- topTable(efit, number = Inf)
  
  #for input
  res <- list()
  res[["DecisionTable"]] <- dt
  res[["DEgenes"]] <- DEgenes
  res[["LogCPM"]] <- v$E
  
  return(res)
}
###############

###############
#function for gene set enrichment analysis (GSEA); TF - do for TFs or no?
#terms is required - it's a gmt file with TF targets or Pathway terms
enrichmentAnalysis <- function(DEgenes, TF = F, terms = NULL, minsize){
  
  gseaInputGeneScores <- DEgenes %>%
    dplyr::select("logFC", "adj.P.Val", "gene")
  
  gseaInputGeneScores$Sign <- sapply(gseaInputGeneScores$logFC, sign)
  gseaInputGeneScores$LogPvalue <- sapply(gseaInputGeneScores$adj.P.Val, function(x){ -log10(x) })
  gseaInputGeneScores$Score <- gseaInputGeneScores$Sign * gseaInputGeneScores$LogPvalue 
  
  #create ranked list of genes
  #gseaInputGeneScores <- DEgenes %>% 
    #mutate(absolute_logFC = abs(logFC)) %>% 
    #dplyr::select(gene, t) %>% 
    #na.omit() %>% 
    #as.data.frame()
  
  genes <- gseaInputGeneScores$gene
  gseaInputGeneScores <- gseaInputGeneScores$Score
  names(gseaInputGeneScores) <- genes
  gseaInputGeneScores <- sort(gseaInputGeneScores, decreasing = T)
  #write(genes, file="~/Desktop/genes.txt")
  
  #perform GSEA 
  fgseaRes <- fgsea(terms, gseaInputGeneScores, minSize=minsize, maxSize=300, nperm=10000)
  
  #filter only statistically significant terms and look at direction (NES > 0 means enrichment at the top of ranked list)
  sign_pos <- fgseaRes %>% filter(padj < 0.05 & NES > 0)
  sign_neg <- fgseaRes %>% filter(padj < 0.05 & NES < 0)
  
  #sort by NES (normalized enrichment score)
  sign_pos <- sign_pos %>% arrange(desc(NES))
  sign_neg <- sign_neg %>% arrange(NES)
  
  #saving resulting data frames to list
  res <- list()
  res[["Genes"]] <- genes
  res[["Enrichment"]] <- fgseaRes
  res[["Scores"]] <- gseaInputGeneScores
  res[["Sign_Pos"]] <- sign_pos
  res[["Sign_Neg"]] <- sign_neg
  return(res)
}
###############


#######################################################
#MAIN DGE FUNCTION
#######################################################
performDGEanalysis <- function(expr_mat, meta_data, group1, group2, termGO = NULL,
                               termTF = NULL, termH = NULL){
  #perform DGE analysis with limma
  test <- limma_DGE_group_analysis(expr_mat, meta_data, group1, group2)
  print("Done with limma analysis")
  #get decision table
  dt <- test[["DecisionTable"]]
  #get data frame with gene expression stats
  DEgenes <- test[["DEgenes"]]
  #get log cpm matrix
  log_cpm <- test[["LogCPM"]]
  
  #arrange genes by logFC
  DEgenes <- DEgenes %>% rownames_to_column("gene") %>% arrange(desc(logFC))
  
  #df of significantly up regulated genes
  upRegulatedGenes <- DEgenes %>% filter(logFC > 1 & adj.P.Val < 0.05)
  #df of significantly down regulated genes
  downRegulatedGenes <- DEgenes %>% filter(logFC < -1 & adj.P.Val < 0.05)
  
  enrichedGoPathways <- enrichmentAnalysis(DEgenes, TF = F, terms = termGO, minsize = 15)
  print("Done with Pathways enrichment")
  sign_pos_pathways <- enrichedGoPathways[["Sign_Pos"]]
  sign_neg_pathways <- enrichedGoPathways[["Sign_Neg"]]
  
  enrichedTFs <- enrichmentAnalysis(DEgenes, TF = T, terms = termTF, minsize = 5)
  print("Done with TFs enrichment")
  sign_pos_TF <- enrichedTFs[["Sign_Pos"]]
  sign_neg_TF <- enrichedTFs[["Sign_Neg"]]
  TFfgsea <- enrichedTFs[["Enrichment"]]
  TFscores <- enrichedTFs[["Scores"]]
  
  enrichedHPathways <- enrichmentAnalysis(DEgenes, TF = F, terms = termH, minsize = 15)
  print("Done with Hallmark Pathways enrichment")
  sign_pos_pathwaysH <- enrichedHPathways[["Sign_Pos"]]
  sign_neg_pathwaysH <- enrichedHPathways[["Sign_Neg"]]
  
  #genes for GOrilla analysis (or any other GO enrichment tool)
  genes <- enrichedGoPathways[["Genes"]]
  
  #saving results
  res <- list()
  res[["DEgenesTable"]] <- DEgenes
  res[["upRegulatedGenes"]] <- upRegulatedGenes
  res[["downRegulatedGenes"]] <- downRegulatedGenes
  res[["Sign_Pos_Pathway"]] <- sign_pos_pathways
  res[["Sign_Neg_Pathway"]] <- sign_neg_pathways
  res[["Sign_HALLMARK_Pos_Pathway"]] <- sign_pos_pathwaysH
  res[["Sign_HALLMARK_Neg_Pathway"]] <- sign_neg_pathwaysH
  res[["Sign_Pos_TF"]] <- sign_pos_TF
  res[["Sign_Neg_TF"]] <- sign_neg_TF
  res[["TF_enrichment"]] <- TFfgsea
  res[["TF_scores"]] <- TFscores
  res[["genes"]] <- genes
  res[["DecisionTable"]] <- dt
  res[["LogCPM"]] <- log_cpm
  return(res)
}

#############################################################################################################
#GET THE MOST VARIABLE GENES
#############################################################################################################
library(DESeq); library(statmod); library(pcaMethods); library(fastICA)

#To get normalized expression magnitudes
lib.size <- estimateSizeFactorsForMatrix(x$counts)
ed <- t(t(x$counts)/lib.size)

#Calculate estimates of variance, coefficient of variation
means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2

#Now fit a regression line based on the controls:
minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
df <- ncol(ed) - 1

#Rank genes by the significance of deviation from the fit
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- ed[varorder,]

#We can also evaluate statistical significance of the deviation
pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
#table(sigVariedGenes)

#heatmap construction
batch <- meta_data$Donor
log_cpm <- cpm(x, log=TRUE)
lcpm.no.batches <- removeBatchEffect(log_cpm, batch)

log_cpm_f <- lcpm.no.batches
colnames(log_cpm_f) <- meta_data$SampleName

#NGFR_unstim, PD1, 41BB, TNFR2, NGFR, OX40, GITR, C4mut, 3zeta, C4wt, NGFR_Tconv, CD28wt_Tconv, NGFR_unstim_Tconv, ICOS, CD28mut, CD28wt
col.order <- c("9649.NGFR_unstim", "960.NGFR_unstim", "967.NGFR_unstim", "9649.PD1", "960.PD1", "967.PD1", "9649.41BB", "960.41BB", "967.41BB",
               "9649.TNFR2", "960.TNFR2", "967.TNFR2", "9649.OX40", "960.OX40", "967.OX40", 
               "9649.GITR", "960.GITR", "967.GITR", "9649.C4mut", "960.C4mut", "967.C4mut", "9649.3zeta", "960.3zeta", "967.3zeta",
               "9649.C4wt", "960.C4wt", "967.C4wt",
               "9649.ICOS", "960.ICOS", "967.ICOS", "9649.CD28mut", "960.CD28mut", "967.CD28mut", 
               "9649.CD28wt", "960.CD28wt", "967.CD28wt", "9649.NGFR", "960.NGFR", "967.NGFR")

samples <- c("NGFR_unstim", "NGFR_unstim", "NGFR_unstim", "PD1", "PD1", "PD1", "41BB", "41BB", "41BB",
               "TNFR2", "TNFR2", "TNFR2", "OX40", "OX40", "OX40", 
               "GITR", "GITR", "GITR", "C4mut", "C4mut", "C4mut", "3zeta", "3zeta", "3zeta",
               "C4wt", "C4wt", "C4wt", 
               "ICOS", "ICOS", "ICOS", "CD28mut", "CD28mut", "CD28mut", "CD28wt", "CD28wt", "CD28wt",
               "NGFR", "NGFR", "NGFR")

log_cpm_f <- log_cpm_f[,col.order]

varGenes <- log_cpm_f[rownames(oed[1:25,]),]
#scaling the data (without it in pheatmap specify scale = "row")
varGenes <- t(scale(t(varGenes)))

my_sample_col <- data.frame(sample = samples)
row.names(my_sample_col) <- colnames(log_cpm_f)

#creating a heatmap with pheatmap (no clustering because we ordered samples by time point)
devSVG(file = "~/Nick_Treg_project/SVG_heatmaps/heatmap_RB_var_25ann.svg")
pheatmap(varGenes, cluster_rows = T, cluster_cols = F, 
         color=bluered(21), scale = "none", clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", show_colnames = T, show_rownames = T, 
         annotation_col = my_sample_col,
         main = "25 most variable genes")
dev.off()

#############################################################################################################
#EXAMPLE OF THE ANALYSIS
#############################################################################################################

#Let's say we want to answer - Which transcriptional pathways are differentially activated in CARs that work in vivo vs. those that don’t.
#Our groups - (CD28wt, CD28mut) AND (3zeta, C4mut, C4wt, GITR, OX40, PD1)

#1A
#1B
#2A
#3A
#4A
#5A
#5B
#6A
groupAll <- c("PD1", "41BB", "TNFR2", "NGFR", "OX40", "GITR", "C4mut", "3zeta", "C4wt", "ICOS", "CD28mut", "CD28wt")
groupCtrl <- c("NGFR_unstim") 
DGEres <- performDGEanalysis(x, meta_data, groupAll, groupCtrl, termGO = termGO,
                             termTF = termTF, termH = termH)

#table with all genes and their scores
DEgenes <- DGEres[["DEgenesTable"]]
#up regulated genes in group1A1 compared to group1A2
upRegulatedGenes <- DGEres[["upRegulatedGenes"]]
#down regulated genes in group1A1 compared to group1A2
downRegulatedGenes <- DGEres[["downRegulatedGenes"]]
#pathways positively enriched in group1A1 compared to group1A2
sign_pos_pathways <- DGEres[["Sign_Pos_Pathway"]]
#pathways negatively enriched in group1A1 compared to group1A2
sign_neg_pathways <- DGEres[["Sign_Neg_Pathway"]]

#pathways positively enriched in group1A1 compared to group1A2
sign_pos_pathwaysH <- DGEres[["Sign_HALLMARK_Pos_Pathway"]]
#pathways negatively enriched in group1A1 compared to group1A2
sign_neg_pathwaysH <- DGEres[["Sign_HALLMARK_Neg_Pathway"]]

TFgsea <- DGEres[["TF_enrichment"]]
TFscores <- DGEres[["TF_scores"]]

#Transcription factors positively enriched in group1A1 compared to group1A2
sign_pos_TF <- DGEres[["Sign_Pos_TF"]]
#Transcription factors negatively enriched in group1A1 compared to group1A2
sign_neg_TF <- DGEres[["Sign_Neg_TF"]]
#Ranked list of genes for GOrilla analysis
genes <- DGEres[["genes"]]
write(genes, file="~/Desktop/genes.txt") #save them to file and provide as input for GOrilla service (http://cbl-gorilla.cs.technion.ac.il/)
#Decision table (potential use - Venn Diagramm; checking what genes are DE in different conditions)
dt <- DGEres[["DecisionTable"]]
#log cpm matrix for heat map construction
log_cpm <- DGEres[["LogCPM"]]

#saving data
write.table(upRegulatedGenes, file="~/Nick_Treg_project/DE_genes/upRegulated_genes_group6A.txt", sep="\t", row.names = F)
write.table(downRegulatedGenes, file="~/Nick_Treg_project/DE_genes/downRegulated_genes_group6A.txt", sep="\t", row.names = F)

write.table(sign_pos_TF[,-8], file="~/Nick_Treg_project/TF/posTF_group6A.txt", sep="\t", row.names = F)
write.table(sign_neg_TF[,-8], file="~/Nick_Treg_project/TF/negTF_group6A.txt", sep="\t", row.names = F)

write.table(sign_pos_pathwaysH[,-8], file="~/Nick_Treg_project/PathwaysHallmark/posPathways_group6A.txt", sep="\t", row.names = F)
write.table(sign_neg_pathwaysH[,-8], file="~/Nick_Treg_project/PathwaysHallmark/negPathways_group6A.txt", sep="\t", row.names = F)

write.table(sign_pos_pathways[,-8], file="~/Nick_Treg_project/Pathways/posPathways_group6A.txt", sep="\t", row.names = F)
write.table(sign_neg_pathways[,-8], file="~/Nick_Treg_project/Pathways/negPathways_group6A.txt", sep="\t", row.names = F)

################################################################################
#heatmap construction
#filtering meta data to groups only
meta_data_f <- meta_data %>% filter(Construct %in% c(groupAll, groupCtrl))
#topGenesRegulated <- cleaned_log_cpm_df[common_regulated_genes,]

#correct for batch - in this case it's donor factor
#batch <- x$samples$Donor
batch <- meta_data_f$Donor
lcpm.no.batches <- removeBatchEffect(log_cpm, batch)

log_cpm_f <- lcpm.no.batches
colnames(log_cpm_f) <- meta_data_f$SampleName

upRegulatedGenes <- upRegulatedGenes$gene
downRegulatedGenes <- downRegulatedGenes$gene
common_regulated_genes <- c(upRegulatedGenes[1:25], downRegulatedGenes[1:25])
topGenesRegulated <- log_cpm_f[common_regulated_genes,]
topGenesRegulated <- topGenesRegulated[,col.order]
#scaling the data (without it in pheatmap specify scale = "row")
topGenesRegulated <- t(scale(t(topGenesRegulated)))

#creating a heatmap with pheatmap (no clustering because we ordered samples by time point)
devSVG(file = "~/Nick_Treg_project/SVG_heatmaps/heatmap_RB_DE_all_ctrl_50.svg")
pheatmap(topGenesRegulated, cluster_rows = T, cluster_cols = F, 
         color=bluered(21), scale = "none", clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", show_colnames = T, show_rownames = T, 
         annotation_col = my_sample_col,
        main = "Clustering heatmap for 50 common regulated genes")
dev.off()

devSVG(file = "~/Nick_Treg_project/SVG_heatmaps/heatmap_6A.svg")
pheatmap(topGenesRegulated, cluster_rows = T, cluster_cols = T, 
         scale = "none", clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", show_colnames = T, show_rownames = FALSE, 
         main = "Clustering heatmap for top common regulated genes")
dev.off()

################################################################################
#Volcano plot construction
devSVG(file = "~/Nick_Treg_project/SVG_volcano/volcano_6A.svg")
ggplot(data = DEgenes, aes(x = logFC, y = -log(adj.P.Val), color = ((-log(adj.P.Val) > 3) & (logFC > 1 | logFC < -1))))+
  scale_colour_manual(name = 'BH p-value < 0.05', values = setNames(c('red','black'),c(T, F)), labels = c("False", "True"))+
  geom_point()+
  geom_vline(xintercept=0)+
  geom_vline(xintercept=-1)+
  geom_vline(xintercept=1)+
  geom_hline(yintercept=3)+
  #geom_text(aes(label = ifelse((-log(adj.P.Val) > 3) & (logFC > 1 | logFC < -1), gene, "")), vjust=-1, size = 3)+
  geom_text(aes(label = ifelse((-log(adj.P.Val) > 5) | (logFC > 2.5 | logFC < -2.5), gene, "")), vjust=-1, size = 3)+#(to visualize some of the genes, note overlaps)
  #geom_text(aes(label = ifelse((-log(adj.P.Val) > 10) & (logFC > 3 | logFC < -3), gene, "")), vjust=-1, size = 3)+
  #xlim(-1.5,1.5)+
  ylab("-log(p-value)")+
  xlab("logFC")+
  labs(title="Gene expression differences in two groups")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size=14),
        strip.background = element_rect(colour="white", fill="white"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))
dev.off()
