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
library(dplyr)
library(RColorBrewer)
library(tibble)
library(ermineR)
library(pheatmap)
library(sva)
library(Rtsne)
library(fgsea)

#Loading the data (change folders to apropriate paths)
files <- list.files("~/Nick_Treg_project/SigDom RNAseq data/RNASeq counts/")
x <- readDGE(files, path = "~/Nick_Treg_project/SigDom RNAseq data/RNASeq counts/", columns=c(1,2)) 
x$counts <- x$counts[1:(nrow(x$counts)-2),]

#Loading metadata and specifying sex factor variable
meta_data <- read.table(file = "~/Nick_Treg_project/SigDom RNAseq data/CAR-Treg RNAseq sample info.csv", sep=",", header=T, row.names = 1)
meta_data$SampleName <- paste(meta_data$Donor, meta_data$Construct, sep=".")
sex <-meta_data$Donor
sex[sex %in% c(9649, 967)] <- "F"
sex[sex %in% c(960)] <- "M"
meta_data$sex <- sex
group <- meta_data$Construct
x$samples$group <- group

#Filtering lowly expressed genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE] #13818    48

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
#1 Which transcriptional pathways are differentially activated in CARs that work in vivo vs. those that don’t.
# Analysis A: (remove CARs that destabilize Treg phenotype after stim)
group1A1 <- c("CD28wt", "CD28mut") #Treg
group1A2 <- c("3zeta", "C4mut", "C4wt", "GITR", "OX40", "PD1") #Treg
# Analysis B: (all constructs)
group1B1 <- c("CD28wt", "CD28mut") #Treg
group1B2 <- c("3zeta", "41BB", "C4mut", "C4wt", "GITR", "ICOS", "OX40", "PD1", "TNFR2") #Treg

#2 Why do 41BB/TNFR2 destabilize Treg phenotype after CAR stimulation?
#Analysis A: (top performing CARs)
group2A1 <- c("41BB", "TNFR2") #Treg
group2A2 <- c("CD28wt", "CD28mut") #Treg
#Analysis B: (all constructs)
group2B1 <- c("41BB", "TNFR2") #Treg
group2B2 <- c("3zeta", "C4mut", "C4wt", "GITR", "ICOS", "OX40", "PD1", "CD28wt", "CD28mut") #Treg

#3 Which transcriptional pathways are differentially activated in a TCR stim vs a CAR stim
#Analysis A: (Tregs)
group3A1 <- c("CD28wt") #Treg
group3A2 <- c("NGFR") #Treg
#Analysis B: (Tconv)
group3B1 <- c("CD28wt_Tconv") #Tconv
group3B2 <- c("NGFR_Tconv") #Tconv

#4 Which transcriptional pathways are differentially activated in Treg vs Tconv CAR stimulation?
group4A1 <- c("CD28wt") #Tconv
group4A2 <- c("CD28_Tconv") #Tconv

#Loading gmt files for enrichment analysis
termGO <- gmtPathways("/home/german/Nick_Treg_project/Human_GOBP_AllPathways_no_GO_iea_April_01_2019_symbol.gmt")
termTF <- gmtPathways("/home/german/Nick_Treg_project/RegNetworkStrong.gmt")

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
  
  #create ranked list of genes
  gseaInputGeneScores <- DEgenes %>% 
    #mutate(absolute_logFC = abs(logFC)) %>% 
    dplyr::select(gene, logFC) %>% 
    na.omit() %>% 
    as.data.frame()
  
  genes <- gseaInputGeneScores$gene
  gseaInputGeneScores <- gseaInputGeneScores$logFC
  names(gseaInputGeneScores) <- genes
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
  res[["Sign_Pos"]] <- sign_pos
  res[["Sign_Neg"]] <- sign_neg
  return(res)
}
###############


#######################################################
#MAIN DGE FUNCTION
#######################################################
performDGEanalysis <- function(expr_mat, meta_data, group1, group2, termGO = NULL,
                               termTF = NULL){
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
  
  #genes for GOrilla analysis (or any other GO enrichment tool)
  genes <- enrichedGoPathways[["Genes"]]
  
  #saving results
  res <- list()
  res[["upRegulatedGenes"]] <- upRegulatedGenes
  res[["downRegulatedGenes"]] <- downRegulatedGenes
  res[["Sign_Pos_Pathway"]] <- sign_pos_pathways
  res[["Sign_Neg_Pathway"]] <- sign_neg_pathways
  res[["Sign_Pos_TF"]] <- sign_pos_TF
  res[["Sign_Neg_TF"]] <- sign_neg_TF
  res[["genes"]] <- genes
  res[["DecisionTable"]] <- dt
  res[["LogCPM"]] <- log_cpm
  return(res)
}

#############################################################################################################
#EXAMPLE OF THE ANALYSIS
#############################################################################################################

#Let's say we want to answer - Which transcriptional pathways are differentially activated in CARs that work in vivo vs. those that don’t.
#Our groups - (CD28wt, CD28mut) AND (3zeta, C4mut, C4wt, GITR, OX40, PD1)

group1A1 <- c("CD28wt", "CD28mut") #Treg
group1A2 <- c("3zeta", "C4mut", "C4wt", "GITR", "OX40", "PD1") #Treg

DGEres <- performDGEanalysis(x, meta_data, group1A1, group1A2, termGO = termGO,
                             termTF = termTF)

#up regulated genes in group1A1 compared to group1A2
upRegulatedGenes <- DGEres[["upRegulatedGenes"]]
#down regulated genes in group1A1 compared to group1A2
downRegulatedGenes <- DGEres[["downRegulatedGenes"]]
#pathways positively enriched in group1A1 compared to group1A2
sign_pos_pathways <- DGEres[["Sign_Pos_Pathway"]]
#pathways negatively enriched in group1A1 compared to group1A2
sign_neg_pathways <- DGEres[["Sign_Neg_Pathway"]]
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

################################################################################
#heatmap construction
#filtering meta data to groups only
meta_data_f <- meta_data %>% filter(Construct %in% c(group1A1, group1A2))
#topGenesRegulated <- cleaned_log_cpm_df[common_regulated_genes,]

#correct for batch - in this case it's donor factor
#batch <- x$samples$Donor
batch <- meta_data_f$Donor
lcpm.no.batches <- removeBatchEffect(log_cpm, batch)

log_cpm_f <- lcpm.no.batches
colnames(log_cpm_f) <- meta_data_f$SampleName

upRegulatedGenes <- upRegulatedGenes$gene
downRegulatedGenes <- downRegulatedGenes$gene
common_regulated_genes <- c(upRegulatedGenes, downRegulatedGenes)
topGenesRegulated <- log_cpm_f[common_regulated_genes,]
#scaling the data (without it in pheatmap specify scale = "row")
topGenesRegulated <- t(scale(t(topGenesRegulated)))

#creating a heatmap with pheatmap (no clustering because we ordered samples by time point)
pheatmap(topGenesRegulated, cluster_rows = T, cluster_cols = T, scale = "none", clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", show_colnames = T, show_rownames = FALSE, 
         main = "Clustering heatmap for top common regulated genes")