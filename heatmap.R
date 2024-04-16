library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(tidyr)
library(edgeR)
library(ComplexHeatmap)
library(stringr)
library("circlize")

###
#Metadata

metadata <- as.data.frame(read_excel("Plin2sh_Ntsh_metafile.xlsx"))
metadata[1:5,1:4]
dim(metadata) 
colnames(metadata)

# combine the Type and Infection columns
metadata$Type_infec=paste(metadata$Type,metadata$Infection,sep = "_")

table(metadata$Infection) 
table(metadata$Batch)
table(metadata$Type_infec)
####


# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Feature_counts_matrix.csv"))
count[1:5,1:5]
count=count[,-1] # remove extra s.no. column
count[1:5,1:5]
count <- tibble::column_to_rownames(count, var = "Geneid")
dim(count) 

length(intersect(colnames(count),metadata$Samples)) 
setdiff(metadata$Samples,colnames(count))  

#######

# Now calculate cpm  
library(edgeR)
cpm <- cpm(count)
cpm[1:5,1:5]
#SMALLEST POSSIBLE GROUP 3
is.exprs <- rowSums(cpm>1) >= 3 # Apply filtering criteria

counts2 <- count[is.exprs, ]
dim(counts2) 

# make samples same in both the files 
rownames(metadata)<- metadata$Samples
common= intersect(rownames(metadata), colnames(counts2))
length(common)  
metadata[1:5,1:5]
counts2[1:5,1:5]
com_exp_mat=counts2[,common]
com_metafile=metadata[common,]

dim(com_exp_mat) 
dim(com_metafile)

table(com_metafile$Batch)
table(com_metafile$Infection)

x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)

vMat <- v$E
dim(vMat) 
vMat[1:5,1:5] 
# In variance_stabilised_counts the rows are the variables and the columns correspond to the samples

# Adjust batch by sva combat here
com_metafile$Batch=as.factor(com_metafile$Batch)
exprMat <- sva::ComBat(vMat, batch=com_metafile$Batch)
exprMat[1:5,1:5]
dim(exprMat)


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
#######################

u1 <- as.data.frame(read_excel("UP_NTsh_unique.xlsx"))
u2 <- as.data.frame(read_excel("UP_Plin2sh_unique.xlsx"))
d1 <- as.data.frame(read_excel("DN_NTsh_unique.xlsx"))
d2 <- as.data.frame(read_excel("DN_Plin2sh_unique.xlsx"))
c1 <- as.data.frame(read_excel("Common_DN_NTsh_Plin2sh.xlsx"))
c2 <- as.data.frame(read_excel("Common_UP_NTsh_Plin2sh.xlsx"))

u1m=merge(exprMat_with_genes_nodup,u1,by="Genes")
u2m=merge(exprMat_with_genes_nodup,u2,by="Genes")
d1m=merge(exprMat_with_genes_nodup,d1,by="Genes")
d2m=merge(exprMat_with_genes_nodup,d2,by="Genes")
c1m=merge(exprMat_with_genes_nodup,c1,by="Genes")
c2m=merge(exprMat_with_genes_nodup,c2,by="Genes")

rownames(u1m)=u1m[,"Genes"]
rownames(u2m)=u2m[,"Genes"]
rownames(d1m)=d1m[,"Genes"]
rownames(d2m)=d2m[,"Genes"]
rownames(c1m)=c1m[,"Genes"]
rownames(c2m)=c2m[,"Genes"]

u1mm = as.matrix(t(u1m[c(-1,-2)]))
u2mm = as.matrix(t(u2m[c(-1,-2)]))
d1mm = as.matrix(t(d1m[c(-1,-2)]))
d2mm = as.matrix(t(d2m[c(-1,-2)]))
c1mm = as.matrix(t(c1m[c(-1,-2)]))
c2mm = as.matrix(t(c2m[c(-1,-2)]))


################ HEATMAP for all g1, g2, g3, g4. g5 and g6 genes ################################################
cbnd_map = cbind(u1mm, u2mm, d1mm, d2mm, c1mm, c2mm)
b= scale(cbnd_map, center=TRUE, scale = TRUE)
cbnd_map[1:4,1:4]
write.csv(cbnd_map, "heatmap_matrix.csv")

cell_colors <- c("NTsh"= "lawngreen","Plin2sh"= "royalblue3")
inf_colors <- c("Inf"="lightpink","Un"="tomato3")
batch_colors <- c("E2"= "lemonchiffon4","E4"= "sandybrown", "E5"= "mediumseagreen")
col_fun = colorRamp2(c(-3,-2, -1,0, 1, 2,3), c("darkblue","blue3","blue1" , "white","red1","red3", "darkred"))

com_metafile$Type=as.factor(com_metafile$Type)
com_metafile$Infection=as.factor(com_metafile$Infection)
com_metafile$Batch=as.factor(com_metafile$Batch)

h1=Heatmap(com_metafile$Type,cluster_rows = F,width = unit(2, "mm"),name="Cell Type",show_row_names = F,col=cell_colors)
h2=Heatmap(com_metafile$Infection,cluster_rows = F,width = unit(2, "mm"),name="Infection status",show_row_names = F,col=inf_colors)
h3=Heatmap(com_metafile$Batch,cluster_rows = F,width = unit(2, "mm"),name="Batch",show_row_names = F,col=batch_colors)
h4=Heatmap(b,name="Scaled Expression",column_names_gp = grid::gpar(fontsize = 3),row_names_gp = grid::gpar(fontsize = 5),col = col_fun)

h1+h2+h3+h4

############################### HEATMAP for all g1, g2, g3, g4. g5 and g6 genes ##################

cbnd_map_2 = cbind(u1mm, u2mm, d1mm, d2mm)
b2= scale(cbnd_map_2, center=TRUE, scale = TRUE)
cbnd_map_2[1:4,1:4]

h1=Heatmap(com_metafile$Type,cluster_rows = F,width = unit(2, "mm"),name="Cell Type",show_row_names = F,col=cell_colors)
h2=Heatmap(com_metafile$Infection,cluster_rows = F,width = unit(2, "mm"),name="Infection status",show_row_names = F,col=inf_colors)
h3=Heatmap(com_metafile$Batch,cluster_rows = F,width = unit(2, "mm"),name="Batch",show_row_names = F,col=batch_colors)
h5=Heatmap(b2,name="Scaled Expression",column_names_gp = grid::gpar(fontsize = 3),row_names_gp = grid::gpar(fontsize = 5),col = col_fun)

h1+h2+h3+h5


###### end ####

