## RNA-seq Differential expression analysis and visualization using DESeq2

# created: 20180326 by Katie (from script written for analysis of Romy's data)
# edited: 20190423 by Nick 

## Input: raw count data
## Output: List of differentially expressed genes for pathway analysis, PCA & heatmap analysis of samples

# sources used: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
# http://www.bioconductor.org/help/workflows/rnaseqGene/

########################################################################################################*
######################################## HOW TO USE ####################################################
########################################################################################################*
##### * Set working directory files in sections 1a and 1b, run sections 0-4a.
##### * After section 4a, the data can be subsetted and exploratory data anlysis plots can be generated in section 5.
##### * Most of the packages listed in section 0 can be found on Bioconductor
########################################################################################################*
########################################################################################################*
########################################################################################################*

####### 0. Reset work envirnment and load required packages #######

# Clear environment
rm(list=ls())

#### Load the libraries required for processing ####
library("DESeq2")
library("biomaRt")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("ggplot2")
library("genefilter")
library("ggrepel")
library("dplyr")
library("data.table")
library("calibrate") 
library(sva)
library(RSvgDevice)



####### 1a. Read raw count data from .txt file #######

# Set working directory to location of files
setwd("~/Desktop/20190305_Nick's RNAseq data/RNASeq counts")

# Create list of text files
txt_files_list = list.files(path="~/Desktop/20190305_Nick's RNAseq data/RNASeq counts", pattern="*.txt") 

# Read the files in, assuming space separator
txt_files_df <- lapply(txt_files_list, function(x) {read.table(file = x, header = F, row.names = 1, sep = "")})

# Combine them
combined_df <- do.call("cbind", lapply(txt_files_df, as.data.frame)) 
colnames(combined_df) = txt_files_list
write.csv(as.data.frame(combined_df), file="CAR-Treg RNAseq raw data.csv") #Write results to csv file
raw_counts = combined_df

# Clean up: Remove last two rows of raw_counts (rows with no feature and ambiguous)
n = nrow(raw_counts)
raw_counts = raw_counts[1:(n-2),]

####### 1b. Import sample property matrix  #######
setwd("~/Desktop/20190305_Nick's RNAseq data/")

# import sample attributes
sampleinfo <- read.csv("CAR-Treg RNAseq sample info.csv", header = TRUE, row.names = 1)

# Make vectors of sample names and corresponding conditions 
donor <- factor(sampleinfo$Donor)
stimulus <- factor(sampleinfo$Stimulus)
cell <- factor(sampleinfo$Cell)
construct <- factor(sampleinfo$Construct)
batch <- factor(sampleinfo$Batch)

# Make a property table of sample name and conditions
SeqDesign <- data.frame(row.names = colnames(raw_counts), donor, construct, cell, batch)

####### 2. Make DESeq data object #######

# Use the DESeqDataSetFromMatrix() function to make a new DESeq data set object
dds <- DESeqDataSetFromMatrix (countData = raw_counts, colData = SeqDesign, design = ~  donor + construct) 
nrow (dds)

# Filter to remove rows with 0 (no reads across all conditions)
# Filter of at least 9 reads is used to try to find genes with at least 1 read for each sample
dds <- dds[ rowSums(counts(dds)) >48, ] 
nrow (dds)

##### 3. Run DESeq, normalize data using rlog transform #####

dds <- DESeq (dds)
dds_rlog <- rlog(dds)
allcounts <- as.matrix(assay(dds_rlog))

##### 4. Apply batch correction and subset data ##### 
# How to use this section: run the section 4a of the code below to prepare data for further exploratory analysis.
# To keep all samples in the analysis, run section 4a, then skip to section 5.
# For subsetting the samples, a filter is first applied to the sample property data frame "SeqDesign", then applied to the combat matrix (which contains all the batch-corrected gene exp data)
# The "combat" object is ultimately the object that is reset and used to generate downstream plots, etc.

#### 4a. Normalize data, extract data and define functions to reset matrices when new subsets are wanted ####
# Applies a 'regularized log' transformation --> transforms count data to log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size

# define reset function
resetmodelmatrix <- function() {
  SeqDesign$construct <<- factor(SeqDesign$construct)
  mod <<- model.matrix(~ as.factor(construct), data = SeqDesign)  # set up design dataset
  mod0 <<- model.matrix(~ 1, data = SeqDesign) # set up design dataset
  print('reset model matrix')
}
resetgeneset <- function() {
  SeqDesign <<- data.frame(row.names = colnames(raw_counts), donor, construct, cell, batch) # resets SeqDesign matrix
  SeqDesign$construct <<- factor(SeqDesign$construct)
  resetmodelmatrix()
  print('reset geneset')
  combat <<- ComBat(allcounts, batch = SeqDesign$batch, mod=mod)
}

resetgeneset() # call batch correction function

#### 4b. Only Tregs ####
resetgeneset() # reset geneset
SeqDesign <- subset(SeqDesign, cell == 'Treg') # select only Tregs
combat <- combat[,colnames(combat) %in% rownames(SeqDesign)]

#### 4c. Only Treg: 41bb, TNFR2, CD28wt/mut  ####
resetgeneset() # reset geneset
constructs <- factor(c('41BB', 'TNFR2', 'CD28wt', 'CD28mut')) # select only constructs of interest
SeqDesign <- subset(SeqDesign, construct %in% constructs) # set up design dataset
combat <- combat[,colnames(combat) %in% rownames(SeqDesign)]

##### 4d. compare to Treg gene signature ######
# this is a list of genes that we have previoulsy published on that delineate Tregs from Tconv and Naive from Memory subsets
resetgeneset()
# SeqDesign <- subset(SeqDesign, cell == 'Treg')
# combat <- combat[,colnames(combat) %in% rownames(SeqDesign)]
Treg.sig <- c("RBMS3", "ICA1", "STAM", "IKZF2", "FOXP3", "LRRC32", "METTL7A", "ZNF532", "VAV3", "ECOP", "TRIB1", "CSF2RB", "HPGD", "TMEM23", "C8ORF70", "CTLA4", "ZBTB38", "TNFRSF1B", "IL1R1", "IL1R2", "IL1RN", "TNFRSF9", "PTPRK", "HDGFRP3", "DACT1", "ID2", "IL7R", "LPIN2", "ABCB1", "ANK3", "NELL2")
combat <- combat[rownames(combat) %in% Treg.sig,]

######### 5. Visualizing the data ########
# How to use this section: After running the subset data operations in section 4, run the section below that generates the desired plot
# To print SVG for publication, run the desired section. An image object, "g", will be generated, then run the section 6 below.

######### 5a. PCA ###### 
combat.pca <- prcomp(t(combat),scale.=T)
library(ggfortify)
tiff(file="PCAgg.tiff", units="in", width=5, height=5, res=300)
g <- autoplot(combat.pca, data = SeqDesign, colour = 'construct', shape = 'donor',
              # frame = TRUE, frame.type = 'construct'
)
print(g)
dev.off()

# tSNE
library(caret)
library(Rtsne)

set.seed(9)
combat.tsne <- Rtsne(t(combat), dims = 2, perplexity = 7)
combat.tsne <- as.data.frame(combat.tsne$Y)

tiff(file="tsne.tiff", units="in", width=7, height=5, res=300)
g <- ggplot(combat.tsne, aes(x=V1, y=V2)) + 
  geom_point(aes(size=0.25, color=SeqDesign$donor)) + 
  geom_text(aes(label=SeqDesign$construct, hjust=0.25, vjust=0.25)) +
  guides(colour=guide_legend(override.aes=list(size=4)),
         size=guide_legend()) +
  xlab("tSNE1") + ylab("tSNE2") +
  ggtitle("t-SNE") +
  theme_light(base_size=8) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
  #scale_colour_brewer(palette = "Set2")
print(g)
dev.off()


######### 5b. Heatmap of all differentially expressed genes ############

# Set colour scale (taken from https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r)
paletteLength <- 50
myColour <- colorRampPalette(c("blue4", "white", "red"))(paletteLength)

# Heatmap of sample distances 
# Plot clustered heatmap of samples & write to file
# (1) Calculate distance matrix from rlog-transformed counts
# This is done using transpose, since dist function expects samples to be rows and different dimensions (genes here) to be columns
SampleDists <- dist(t(combat))
SampleDistMatrix <- as.matrix(SampleDists)

# (2) Label rows and columns with sample names
rownames(SampleDistMatrix) <- paste(SeqDesign$donor, SeqDesign$construct, sep="-") #use SeqDesign2 for only Treg; SeqDesign for all
colnames(SampleDistMatrix) <- NULL

# (3) Plot heatmap and save to .png file
#colours <- colorRampPalette(rev(brewer.pal(9, 'Blues')))(255)
tiff(file="Heatmap_samples.tiff")
g <- pheatmap(SampleDistMatrix, clustering_distance_rows=SampleDists, clustering_distance_cols=SampleDists, col = myColour, clustering_method = "ward.D2")
print(g)
dev.off()

# Heatmap-all genes: takes regularized log transformation (rlog), computes difference between rlog of sample and row mean
tiff(file="Heatmap_allgenes.tiff", units="in", width=4, height=6, res=300)
mat2 <- combat
colnames(mat2) <- paste(SeqDesign$donor, SeqDesign$construct, sep="-") 
mat2 <- mat2 - rowMeans(mat2)
myBreaks <- unique(c(seq(min(mat2), -0.001, length = 24), 0, seq(0.001, max(mat2), length = 25))) #Sets number of breaks in colour scale, with 0 as the mid-point
g <- pheatmap(mat2, color = myColour, breaks = myBreaks, cluster_cols = TRUE, cluster_rows = TRUE, 
              show_rownames = TRUE, 
              treeheight_row = 0, clustering_method = "mcquitty")
print(g)
dev.off ()


####### 5c. Heatmap of selected genes ####### 

# generate ordered p-value data table
resetmodelmatrix()
combat.sortp = f.pvalue(combat,mod,mod0) # generates a p value for each differentially-expressed gene in the model matrix
combat.sortp <- data.frame(names(combat.sortp),combat.sortp)
colnames(combat.sortp) <- c('gene','pvalue')
combat.sortp <- combat.sortp[order(combat.sortp$pvalue),]

## Print heatmap for genes with p < 0.01 ##
# Heatmap of genes of less than a specified p value
# adjust p value of interest and subset genelists
genelist <- subset(combat.sortp, pvalue < 0.01)
combat.select <- subset(combat, rownames(combat) %in% genelist$gene)
tiff(file="Heatmap_p01.tiff", units="in", width=10, height=6, res=300)
mat2 <- combat.select
colnames(mat2) <- paste(SeqDesign$donor, SeqDesign$construct, sep="-") 
mat2 <- mat2 - rowMeans(mat2)
myBreaks <- unique(c(seq(min(mat2), -0.001, length = 24), 0, seq(0.001, max(mat2), length = 25))) #Sets number of breaks in colour scale, with 0 as the mid-point
g <- pheatmap(mat2, show_rownames = FALSE, color = myColour, breaks = myBreaks, cluster_cols = TRUE, cluster_rows = TRUE, 
         treeheight_row = 0, clustering_method = "ward.D2")
print(g)
dev.off ()

## Print heatmap of top X genes in p value list ##
genelist <- combat.sortp[1:50,]
combat.select <- subset(combat, rownames(combat) %in% genelist$gene)
tiff(file="Heatmap_top50.tiff", units="in", width=3, height=6, res=300)
mat2 <- combat.select
# colnames(mat2) <- paste(SeqDesign$donor, SeqDesign$construct, sep="-") 
colnames(mat2) <- SeqDesign$construct
mat2 <- mat2 - rowMeans(mat2)
min(mat2)
max(mat2)
myBreaks <- unique(c(seq(min(mat2), -0.001, length = 24), 0, seq(0.001, max(mat2), length = 25))) #Sets number of breaks in colour scale, with 0 as the mid-point
g <- pheatmap(mat2, show_rownames = TRUE, color = myColour, breaks = myBreaks, cluster_cols = TRUE, cluster_rows = TRUE, 
         treeheight_row = 0, fontsize = 6, clustering_method = "mcquitty")
print(g)
dev.off ()

##### 6. print SVG #####
devSVG(file="file.svg", width=3, height=6)
  print(g)
dev.off()

