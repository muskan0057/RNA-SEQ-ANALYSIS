## DEG FOR 4 different groups ##



# load required libraries

getwd()
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(readxl)
library(qtl2)
library(dplyr)
library(tidyverse)
library(data.table)
library(tidyr)
library(edgeR)


# metafile  
metadata <- as.data.frame(read_excel("Plin2sh_Ntsh_metafile.xlsx"))
metadata[1:5,1:4]
dim(metadata) #  12  4
colnames(metadata)


# combine the Type and Infection columns
metadata$Type_infec=paste(metadata$Type,metadata$Infection,sep = "_")

table(metadata$Infection) 
# Inf  Un 
#  6   6 

table(metadata$Batch)
# E2 E4 E5 
# 4  4  4 

table(metadata$Type_infec)
####



# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Feature_counts_matrix.csv"))
count[1:5,1:5]
count=count[,-1] # remove extra s.no. column
count[1:5,1:5]
count <- tibble::column_to_rownames(count, var = "Geneid")
dim(count)  #  63241    12

length(intersect(colnames(count),metadata$Samples)) #  12  (all samples match)
setdiff(metadata$Samples,colnames(count))  # 0

#######

# Now calculate cpm  
library(edgeR)
cpm <- cpm(count)
write.csv(cpm, "CPM.csv")
cpm[1:5,1:5]
#SMALLEST POSSIBLE GROUP 3
is.exprs <- rowSums(cpm>1) >= 3 # Apply filtering criteria

counts2 <- count[is.exprs, ]
dim(counts2) # 19622    12


# make samples same in both the files 
rownames(metadata)<- metadata$Samples
common= intersect(rownames(metadata), colnames(counts2))
length(common)  # 12
metadata[1:5,1:5]
counts2[1:5,1:5]
com_exp_mat=counts2[,common]
com_metafile=metadata[common,]

dim(com_exp_mat) # 19622    12
dim(com_metafile) # 12  5

table(com_metafile$Batch)
table(com_metafile$Infection)



## EdgeR 
d <- DGEList(counts=com_exp_mat, group=com_metafile$Type_infec)  # group here is combined column

# Normalize the data
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)


## convert to factor
batches <- com_metafile$Batch 
table(batches)
# E2 E4 E5 
#  4  4  4
batches=as.factor(batches)

group <- com_metafile$Type_infec
table(group)
# NTsh_Inf     NTsh_Un  Plin2sh_Inf  Plin2sh_Un 
#       3           3           3           3
group=as.factor(group)


# create modelDesign
modelDesign <- model.matrix(~0 + group + batches )
modelDesign1 <- model.matrix(~group + batches )
modelDesign2 <- model.matrix(~batches + group )
head(modelDesign)
dim(modelDesign) # 12  6

table(modelDesign[,1])
table(modelDesign[,2])
table(modelDesign[,3])

# let us now make all the comparisons at this step

contrast_one <- makeContrasts(
  Plin2shVsNTsh_Inf=groupPlin2sh_Inf-groupNTsh_Inf,
  Plin2shVsNTsh_Un=groupPlin2sh_Un-groupNTsh_Un,
  InfVsUn_NTsh=groupNTsh_Inf-groupNTsh_Un, 
  InfVsUn_Plin2sh=groupPlin2sh_Inf-groupPlin2sh_Un,
  levels=modelDesign)

contrast_one
#                    Contrasts
# Levels             Plin2shVsNTsh_Inf Plin2shVsNTsh_Un InfVsUn_NTsh InfVsUn_Plin2sh
# groupNTsh_Inf                   -1                0            1               0
# groupNTsh_Un                     0               -1           -1               0
# groupPlin2sh_Inf                 1                0            0               1
# groupPlin2sh_Un                  0                1            0              -1
# batchesE4                        0                0            0               0
# batchesE5                        0                0            0               0

nrow(modelDesign)
ncol(d$counts)

fit_glm <- glmFit(d,modelDesign)
ensg_hgnc=read.table("ensg_hgnc",header=F)
head(ensg_hgnc)
rownames(ensg_hgnc)=ensg_hgnc$V1
## Apply loop here
for (i in 1:ncol(contrast_one)){
  onevsrest <- glmLRT(fit_glm , contrast = contrast_one[,i])
  tt_onevsrest <- topTags(onevsrest,n=nrow(d))
  dim(tt_onevsrest) # 19622   8
  toptable=tt_onevsrest$table
  head(toptable)
  dim(toptable) # 19622   8
  d=merge(ensg_hgnc,toptable,by=0)
  head(d)
  d=d[,-1]
  colnames(d)[1]="ENSG"
  colnames(d)[2]="HGNC"
  write.csv(d,file=paste("toptable_",colnames(contrast_one)[i],".csv",sep=""))
  sel_toptable=d[(d$logFC >1.00 |d$logFC < -1.00) & d$FDR < 0.05,]
  write.csv(sel_toptable,file=paste("sel_toptable_",colnames(contrast_one)[i],".csv",sep=""))
  
  print(colnames(contrast_one)[i])
}



####### 
## PCA ## 



# Read the counts matrix
count=as.data.frame(read_csv("Feature_counts_matrix.csv"),check.names = F)
count[1:5,1:5]
count=count[,-1]
colnames(count)
count=tibble::column_to_rownames(count,var = "Geneid")
dim(count)   # 63241    12
count[1:5,1:5]


# Read the metafile
metadata <- as.data.frame(read_excel("Plin2sh_Ntsh_metafile.xlsx"))
metadata[1:5,1:4]
dim(metadata) # 12  4
colnames(metadata)  
table(metadata$Batch)
# E2 E4 E5 
#  4  4  4


# match samples in count matrix and metadata 
# Check the common samples
length(intersect(metadata$Samples,colnames(count))) # 12
setdiff(metadata$Samples,colnames(count))  # 0  --> all samples match

count <- count[,metadata$Samples]
dim(count)  # 63241    12
count[1:5,1:5]

rownames(metadata)=metadata[,1]  # making samples column as rownames
metadata[1:5,1:4]
dim(metadata)  # 12  4
all(colnames(count) %in% metadata$Samples)  # TRUE --> samples in count and metadata are same

#############--------------------------

# 1. PCA by Adjusting batch effect
library(sva)

class(count)
count=as.matrix(count)

count[1:5,1:5]
dim(count)  # 63241    12

cpms <- cpm(count)  # calculate cpm
cpms[1:5,1:5]
dim(cpms) # 63241    12

## Apply filter now
keep <- rowSums(cpms > 1) >= 3
counts2=count[keep,]
dim(counts2) # 19622    12
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=F)

vMat <- v$E

# adjust batch here
exprMat <- sva::ComBat(vMat, batch=metadata$Batch)
dim(exprMat)

metadata$infection_type=paste(metadata$Type,metadata$Infection,sep="_")
metadata$infection_type=as.factor(metadata$infection_type)

sel = order(apply(exprMat, 1, var), decreasing=TRUE)[1:1000]
vMat_sel=exprMat[sel,]
dim(vMat) # 19622  12
dim(vMat_sel) # 1000  12
vMat[1:5,1:5]
t_vMat_sel=t(vMat_sel)
t_vMat_sel[1:5,1:5]
#### 

mypca <- function(x, center = TRUE, scale = TRUE){
  # Samples should be in rows
  # Variables (genes) in the columns
  
  # remove constant variables
  constant_val = apply(x,2,'sd')
  x_reduced = x[,constant_val>0]
  
  # perform SVD
  SVD <- svd(scale(x_reduced,center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frams
  loadings <- data.frame(SVD$v)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  rownames(loadings) <- colnames(x_reduced)
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  
  # return result
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}
pca_results <- mypca(t_vMat_sel, 
                     center = TRUE, 
                     scale = TRUE)
scores <- pca_results$scores
dim(scores) # 12  12
scores[1:5,1:5] # 12 PC
library("dplyr")
library("tidyverse")

scores_with_conditions <- 
  scores %>% 
  rownames_to_column("Samples") %>% # to prepare to join on the "sample" column
  left_join(x = .,                 # this means that we are passing the 'scores' dataframe 
            y = metadata,         # this dataframe contains the sample to condition correspondence
            by = "Samples")
head(scores_with_conditions)
dim(scores_with_conditions)  # 12  16
dim(metadata) # 12  4
dim(scores)  # 12  12
colnames(scores_with_conditions)

## variance explained
explained_variance <- pca_results$explained_var %>% pull("exp_var")
colnames(scores_with_conditions)
table(scores_with_conditions$Infection) 
# Inf  Un 
#  6   6 
dim(scores_with_conditions) # 12  16

scores_with_conditions$Batch=as.factor(scores_with_conditions$Batch)
scores_with_conditions$Infection=as.factor(scores_with_conditions$Infection)
table(scores_with_conditions$Infection)
table(scores_with_conditions$Batch)

colnames(scores_with_conditions)

library("ggplot2")
library("colorspace")
################ PCA PLOT
pdf("FC_1000_withadjust_new.pdf", width = 20, height = 8)
ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = infection_type, shape = Batch)) + scale_shape_manual(values = c(19,8,5)) +
  geom_point(size = 5)  +
  scale_fill_manual(values = c("limegreen","yellow3", "violetred3"))+
  #stat_ellipse() +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) + 
  #theme_dose(font.size = 10)+
  ggtitle("PCA plot for Infected vs Uninfected and Batch overlaid with batch correction") +
  labs(shape="Batch", color="Type of sample")+
  scale_color_manual(values = c("limegreen","purple", "violetred3","orange"))
#scale_color_manual(brewer.pal(3,"Set1"))
# scale_color_manual(values = c("limegreen","lightslateblue","lightsalmon4",
#                               "yellow3", "violetred3","turquoise2","tomato2",
#                               "thistle4","springgreen4","lightpink1","orange","maroon")

dev.off()
getwd()
####################################################################################
# WITHOUT BATCH CORRECTION (NO COMBAT_SEQ)

cpms <- cpm(count)
cpms[1:5,1:5]
dim(cpms) #  63241    12
## Apply filter
keep <- rowSums(cpms > 1) >= 3
counts=count[keep,]
dim(counts) #  19622    12
x <- DGEList(counts=counts)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=F)

vMat <- v$E
sel = order(apply(vMat, 1, var), decreasing=TRUE)[1:1000 ]
vMat_sel=vMat[sel,]
dim(vMat) #  19622    12
dim(vMat_sel) # 1000  12
vMat[1:5,1:5]
t_vMat_sel=t(vMat_sel)


mypca <- function(x, center = TRUE, scale = TRUE){
  # Samples should be in rows
  # Variables (genes) in the columns
  
  # remove constant variables
  constant_val = apply(x,2,'sd')
  x_reduced = x[,constant_val>0]
  
  # perform SVD
  SVD <- svd(scale(x_reduced,center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frams
  loadings <- data.frame(SVD$v)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  rownames(loadings) <- colnames(x_reduced)
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  
  # return result
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}
pca_results <- mypca(t_vMat_sel, 
                     center = TRUE, 
                     scale = TRUE)
scores <- pca_results$scores
dim(scores) # 12  12
scores[1:5,1:5] # 12 PC
library("dplyr")
library("tidyverse")
metadata$infection_type=paste(metadata$Type,metadata$Infection,sep="_")
scores_with_conditions <- 
  scores %>% 
  rownames_to_column("Samples") %>% # to prepare to join on the "sample" column
  left_join(x = .,                 # this means that we are passing the 'scores' dataframe 
            y = metadata,         # this dataframe contains the sample to condition correspondence
            by = "Samples")
head(scores_with_conditions)
dim(scores_with_conditions)  #  12  16
dim(metadata)  # 12  4
dim(scores)  #  12  12
colnames(scores_with_conditions)
## variance explained
explained_variance <- pca_results$explained_var %>% pull("exp_var")
colnames(scores_with_conditions)
table(scores_with_conditions$Type) 
# NTsh Plin2sh 
# 6       6 

scores_with_conditions$Batch=as.factor(scores_with_conditions$Batch)
scores_with_conditions$Infection=as.factor(scores_with_conditions$Infection)
table(scores_with_conditions$Infection)
table(scores_with_conditions$Batch)
# E2 E4 E5 
#  4  4  4 

library("ggplot2")
################ PCA PLOT
pdf("FC_1000_withoutadjust_new.pdf", width = 20, height = 8)
ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = infection_type, shape = Batch)) + scale_shape_manual(values = c(19,8,5)) +
  geom_point(size = 5)  +
  scale_fill_manual(values = c("limegreen","purple", "violetred3","orange"))+
  #stat_ellipse() +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) + 
  #theme_dose(font.size = 10)+
  ggtitle("PCA plot for Infected vs Uninfected and Batch overlaid without batch correction") +
  labs(shape="Batch", color="Type of sample")+
  scale_color_manual(values = c("limegreen","purple", "violetred3","orange"))

dev.off()

############################################



####################################
####################################################

# For heatmap 
library(ComplexHeatmap)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(readr)
library("edgeR")

dim(com_exp_mat) # 19622    12
dim(com_metafile) # 12  5

x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)

vMat <- v$E
dim(vMat) # 19622    12
vMat[1:5,1:5] 
# In variance_stabilised_counts the rows are the variables and the columns correspond to the samples

# Adjust batch by sva combat here
com_metafile$Batch=as.factor(com_metafile$Batch)
exprMat <- sva::ComBat(vMat, batch=com_metafile$Batch)
exprMat[1:5,1:5]
dim(exprMat) # 19622    12

############################--------------------------------------------

## ADD GENE NAMES NOW 

# Now add the gene names corresponding to the ensembl ids
gene_names=scan("./ensg_hgnc",what="",sep="\t")
gene_names<- as.data.frame(gene_names)
dim(gene_names)  # 63241     1

# Split name column into ensembl ids and gene names
library(stringr)
gene_names[c('Ensembl_ids', 'Genes')] <- str_split_fixed(gene_names$gene_names, ' ', 2)
head(gene_names)
gene_names=gene_names[,-1]  # removing first column
head(gene_names)

toberemoved= grep("ENS",gene_names$Genes) # remove the ENSEMBL ids present in Genes column
gene_names2=gene_names[-toberemoved,]
dim(gene_names2) # 42636     2



# Let us now add gene names to the exprMat created above
exprMat=tibble::rownames_to_column(as.data.frame(exprMat),"Ensembl_ids")
exprMat[1:5,1:4]
dim(exprMat) #  19622     13

exprMat_with_genes=merge(exprMat,gene_names2,by="Ensembl_ids")  # merge
exprMat_with_genes=exprMat_with_genes %>% relocate(Genes)
dim(exprMat_with_genes)  # 16382    14

exprMat_with_genes[1:5,1:5]
length(unique(exprMat_with_genes$Genes))  # 16321
table(duplicated(exprMat_with_genes$Genes))
# FALSE  TRUE
# 16321    61

rownames(exprMat_with_genes)=exprMat_with_genes[,2]  # making ensembl ids as rownames

##  To remove duplicate genes
e_g=exprMat_with_genes[,1:2] # assigning gene names and ensembl ids column to a new variable

# Sample DataFrame
e_g$Variance <- apply(exprMat_with_genes[,3: dim(exprMat_with_genes)[2]], 1, var) # calculate variance of each gene
dim(e_g) # 16382     3

e_g[duplicated(e_g$Genes),] # check for duplicate genes

# Sample data frame
library("dplyr")
result <- e_g %>%
  group_by(Genes) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

# Print the result
dim(result) # 16321     3
head(result)


# Only keep these genes with high variance in the matrix
exprMat_with_genes_nodup=exprMat_with_genes[result$Ensembl_ids,]
dim(exprMat_with_genes_nodup) # 16321    14
exprMat_with_genes_nodup[1:5,1:5]

# rownames(exprMat_with_genes_nodup)=exprMat_with_genes_nodup[,1]  # making the gene names as rownames here
# exprMat_with_genes_nodup=exprMat_with_genes_nodup[,3:dim(exprMat_with_genes_nodup)[2]]
# dim(exprMat_with_genes_nodup) # 16321    12
# exprMat_with_genes_nodup[1:5,1:5]

############################---------------------------------------------


# Keep only the common up and downregulated genes in the matrix
sel_exprMat_Plin2shVsNTsh_Inf=exprMat_with_genes_nodup[ gg_Plin2shVsNTsh_Inf,]    # for Plin2shVsNTsh_Inf
dim(sel_exprMat_Plin2shVsNTsh_Inf) # 15  14
sel_exprMat_Plin2shVsNTsh_Inf[1:4,1:5]
sel_exprMat_Plin2shVsNTsh_Inf=na.omit(sel_exprMat_Plin2shVsNTsh_Inf)
dim(sel_exprMat_Plin2shVsNTsh_Inf)  # 14  14


sel_exprMat_Plin2shVsNTsh_Un=exprMat_with_genes_nodup[ gg_Plin2shVsNTsh_Un,]     # for Plin2shVsNTsh_Un
dim(sel_exprMat_Plin2shVsNTsh_Un) # 3  14
sel_exprMat_Plin2shVsNTsh_Un[1:4,1:5]
sel_exprMat_Plin2shVsNTsh_Un=na.omit(sel_exprMat_Plin2shVsNTsh_Un)
dim(sel_exprMat_Plin2shVsNTsh_Un)  # 3  14


sel_exprMat_InfVsUn_NTsh=exprMat_with_genes_nodup[ gg_InfVsUn_NTsh,]             # for InfVsUn_NTsh
dim(sel_exprMat_InfVsUn_NTsh) # 412  14
sel_exprMat_InfVsUn_NTsh[1:4,1:5]
sel_exprMat_InfVsUn_NTsh=na.omit(sel_exprMat_InfVsUn_NTsh)
dim(sel_exprMat_InfVsUn_NTsh)  # 377  14


sel_exprMat_InfVsUn_Plin2sh=exprMat_with_genes_nodup[ gg_InfVsUn_Plin2sh,]       # for InfVsUn_Plin2sh
dim(sel_exprMat_InfVsUn_Plin2sh) # 240  14
sel_exprMat_InfVsUn_Plin2sh[1:4,1:5]
sel_exprMat_InfVsUn_Plin2sh=na.omit(sel_exprMat_InfVsUn_Plin2sh)
dim(sel_exprMat_InfVsUn_Plin2sh)  # 223  14


## Read the gene list given by Sheetal ma'am
# find the common with the obtained up/downreg genes and make heatmap for them 

gene_list <- as.data.frame(read_excel("../List for heatmap.xlsx"))
gene_list_sel <- gene_list[,3]  # keeping the 2nd list from the file which is the longest
gene_list_sel=as.data.frame(gene_list_sel)
colnames(gene_list_sel)=gene_list_sel[1,]
gene_list_sel=gene_list_sel[-1,]
length(gene_list_sel)  # 227  

## keep just the common genes from this list in the matrix

#### 1.
sel_exprMat_Plin2shVsNTsh_Inf[1:5,1:4]
rownames(sel_exprMat_Plin2shVsNTsh_Inf)=sel_exprMat_Plin2shVsNTsh_Inf[,1]  # make gene names as rownames
sel_exprMat_Plin2shVsNTsh_Inf=sel_exprMat_Plin2shVsNTsh_Inf[,-c(1:2)] # remove gene names and ensembl id columns
dim(sel_exprMat_Plin2shVsNTsh_Inf)  # 14  12

sel_exprMat_Plin2shVsNTsh_Inf=sel_exprMat_Plin2shVsNTsh_Inf[rownames(sel_exprMat_Plin2shVsNTsh_Inf) %in% gene_list_sel,]
dim(sel_exprMat_Plin2shVsNTsh_Inf)  # 1  12

#### 2.
sel_exprMat_Plin2shVsNTsh_Un[1:5,1:4]
rownames(sel_exprMat_Plin2shVsNTsh_Un)=sel_exprMat_Plin2shVsNTsh_Un[,1]
sel_exprMat_Plin2shVsNTsh_Un=sel_exprMat_Plin2shVsNTsh_Un[,-c(1,2)]
dim(sel_exprMat_Plin2shVsNTsh_Un)  # 3  12

sel_exprMat_Plin2shVsNTsh_Un=sel_exprMat_Plin2shVsNTsh_Un[rownames(sel_exprMat_Plin2shVsNTsh_Un) %in% gene_list_sel,]
dim(sel_exprMat_Plin2shVsNTsh_Un)  # 0  12

#### 3.
sel_exprMat_InfVsUn_NTsh[1:5,1:4]
rownames(sel_exprMat_InfVsUn_NTsh)=sel_exprMat_InfVsUn_NTsh[,1]
sel_exprMat_InfVsUn_NTsh=sel_exprMat_InfVsUn_NTsh[,-c(1,2)]
dim(sel_exprMat_InfVsUn_NTsh)  # 377  12

sel_exprMat_InfVsUn_NTsh=sel_exprMat_InfVsUn_NTsh[rownames(sel_exprMat_InfVsUn_NTsh) %in% gene_list_sel,]
dim(sel_exprMat_InfVsUn_NTsh)  # 106  12

#### 4.
sel_exprMat_InfVsUn_Plin2sh[1:5,1:4]
rownames(sel_exprMat_InfVsUn_Plin2sh)=sel_exprMat_InfVsUn_Plin2sh[,1]
sel_exprMat_InfVsUn_Plin2sh=sel_exprMat_InfVsUn_Plin2sh[,-c(1,2)]
dim(sel_exprMat_InfVsUn_Plin2sh)  # 223  12

sel_exprMat_InfVsUn_Plin2sh=sel_exprMat_InfVsUn_Plin2sh[rownames(sel_exprMat_InfVsUn_Plin2sh) %in% gene_list_sel,]
dim(sel_exprMat_InfVsUn_Plin2sh)  # 98  12



## HEATMAP -----
t_sel=t(sel_exprMat_InfVsUn_Plin2sh)  # transpose so rows are samples and columns correspond to genes
dim(t_sel)
ti=merge(com_metafile,t_sel,by = 0)

#ti_Plin2shVsNTsh_Inf=ti[!ti$Infection=="Un",]  # removing Un samples from here

#ti_Plin2shVsNTsh_Un=ti[!ti$Infection=="Inf",]  # removing  Inf samples from here

ti_InfVsUn_NTsh=ti[!ti$Type=="Plin2sh",]  # removing  Plin2sh samples from here

ti_InfVsUn_Plin2sh=ti[!ti$Type=="NTsh",]  # removing  NTsh samples from here

ti_InfVsUn_Plin2sh=ti_InfVsUn_Plin2sh[order(ti_InfVsUn_Plin2sh$Type),,drop=F]
rownames(ti_InfVsUn_Plin2sh)=ti_InfVsUn_Plin2sh$Row.names
meta1=ti_InfVsUn_Plin2sh[,1:6]
exp=ti_InfVsUn_Plin2sh[,7:dim(ti_InfVsUn_Plin2sh)[2]]
dim(meta1)  # 6  6
dim(exp)  
exp[1:4,1:4]
meta1[1:4,1:4]


library(ComplexHeatmap)
library("circlize")
Heatmap(scale(exp))

pdf("heatmap_InfVsUn_Plin2sh.pdf",width=8,height = 10)

meta1=meta1[order(meta1$Type),]
exp=exp[rownames(meta1),]
head(meta1)
meta11=meta1[,c(3:5)]
head(meta11)
meta11$Batch

cell_colors <- c("NTsh"= "lawngreen","Plin2sh"= "royalblue3")
inf_colors <- c("Inf"="lightpink","Un"="tomato3")
batch_colors <- c("E2"= "lemonchiffon4","E4"= "sandybrown", "E5"= "mediumseagreen")
col_fun = colorRamp2(c(-3,-2, -1,0, 1, 2,3), c("yellow","yellow3","yellow3" , "black","plum","plum3", "magenta2"))


meta11$Type=as.factor(meta11$Type)
meta11$Infection=as.factor(meta11$Infection)
meta11$Batch=as.factor(meta11$Batch)

h1=Heatmap(meta11$Type,cluster_rows = F,width = unit(1, "cm"),name="Cell Type",show_row_names = F,col=cell_colors)
h2=Heatmap(meta11$Infection,cluster_rows = F,width = unit(1, "cm"),name="Infection status",show_row_names = F,col=inf_colors)
h3=Heatmap(meta11$Batch,cluster_rows = F,width = unit(1, "cm"),name="Batch",show_row_names = F,col=batch_colors)
h4=Heatmap(scale(cbnd_map),name="Scaled Expression",column_names_gp = grid::gpar(fontsize = 3),row_names_gp = grid::gpar(fontsize = 5),col = col_fun)

h1+h2+h3+h4
dev.off()
###############

# pdf("non_cluster_heatmap_NTshVsPlin2sh_Un.pdf")
# h4 + h3 + h2 + h1 
# dev.off()
####################

##---------------------------------------------------------------------------------------------------------


