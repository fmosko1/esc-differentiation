library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(biomaRt)
library(patchwork)

lookup = read.csv("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/data/SC_2023/MH-108-sequencing_ID.csv")

# biomart lookup between ensembl and mgi symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")#, host = "uswest.ensembl.org")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart),]

#load data
seq.data.239246 = Read10X(data.dir = "/cellfile/datapublic/ftitztei/Rscripts/TCpackage/data/SC_2023/multi/S242614_multicount/per_sample_outs/239246/count/sample_filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data)
seq.239246 = CreateSeuratObject(counts = seq.data.239246, project = "Sequencing_Data_239246", min.cells = 0, min.features = 0)

#merge SeuratObjects
seq.combined = merge(seq.239237, y = c(seq.239238, seq.239243, seq.239244, seq.239245, seq.239246), add.cell.ids = c("239237","239238","239243","239244","239245","239246"), project = "Sequencing Data")
head(colnames(seq.combined))
table(seq.combined$orig.ident)

seq.combined <- JoinLayers(seq.combined, assay = "RNA")


### Quality Control ###

mitochondrial = c("ENSMUSG00000064348", "ENSMUSG00000064344",
                   "ENSMUSG00000064365", "ENSMUSG00000064366",
                   "ENSMUSG00000064371", "ENSMUSG00000064339",
                   "ENSMUSG00000064350", "ENSMUSG00000064340",
                   "ENSMUSG00000064349", "ENSMUSG00000064338",
                   "ENSMUSG00000064363", "ENSMUSG00000064367",
                   "ENSMUSG00000064361", "ENSMUSG00000064354",
                   "ENSMUSG00000064358", "ENSMUSG00000064352",
                   "ENSMUSG00000064337", "ENSMUSG00000064359",
                   "ENSMUSG00000064346", "ENSMUSG00000064353",
                   "ENSMUSG00000064357", "ENSMUSG00000064370",
                   "ENSMUSG00000064345", "ENSMUSG00000064364",
                   "ENSMUSG00000064351", "ENSMUSG00000064368",
                   "ENSMUSG00000064369", "ENSMUSG00000064336",
                   "ENSMUSG00000064341", "ENSMUSG00000064355",
                   "ENSMUSG00000064347", "ENSMUSG00000064343",
                   "ENSMUSG00000064356", "ENSMUSG00000064342",
                   "ENSMUSG00000064360", "ENSMUSG00000065947",
                   "ENSMUSG00000064372")

mitochondrial.mgi = biomart[which(biomart$ensembl_gene_id %in% mitochondrial), 3 ]

seq.combined[["percent.mt"]] = PercentageFeatureSet(seq.combined,features = mitochondrial.mgi[which(mitochondrial.mgi %in% row.names(seq.combined))])

seq.combined[["nonmt.libsize"]] <- colSums(seq.combined@assays[["RNA"]]@layers[["counts"]][which(!row.names(seq.combined@assays[["RNA"]]@layers[["counts"]])%in%
                                                                                  mitochondrial.mgi),])
sum(tags[which(tags$X%in%colnames(seq.combined)),1]==colnames(seq.combined))
seq.combined[["type"]] <- tags[which(tags$X%in%colnames(seq.combined)),2]

#Visualize QC metrics as violin plot
VlnPlot(seq.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol = 3), pt.size = 0)
seq.combined = seq.combined[which(lookup$ID %in% seq.combined@active.ident), ]
#Make vlnplot better -> rename samples
test <- seq.combined@active.ident
levels(test)[levels(test) == "Sequencing_Data_239246"] <- "N24_WT_2"
levels(test)
seq.combined@active.ident <- test

#sum(active.ident[which(seq.combined@active.ident %in% colnames(seq.combined)),1]==colnames(seq.combined))
#seq.combined[["type"]] <- tags[which(tags$X%in%colnames(seq.combined)),2]

#ggplot QC metrics
gg_qc = cbind(seq.combined[["nFeature_RNA"]],seq.combined[["nCount_RNA"]],
              seq.combined[["percent.mt"]])
ggplot(gg_qc,aes(x=nCount_RNA,y=percent.mt,
                 alpha = 0.05, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_point() +
  theme_bw()

ggplot(gg_qc,aes(x=nFeature_RNA,y=percent.mt,
                 alpha = 0.05, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_point()+
  theme_bw()

ggplot(gg_qc,aes(x=nCount_RNA,y=nFeature_RNA,
                 alpha = 0.05, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_point()+
  theme_bw()

## Fig 5.21 B
ggplot(gg_qc,aes(x=seq.combined@active.ident,y=nFeature_RNA,
                 alpha = 0.2, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_violin()+
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "grey","black"))+
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                             "darkorange3","coral3",
                             "darkred","grey","black"))+
  theme_bw()

## Fig 5.21 A
ggplot(gg_qc,aes(x=seq.combined@active.ident,y=nCount_RNA,
                 alpha = 0.2, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_violin()+
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "grey","black"))+
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                             "darkorange3","coral3",
                             "darkred","grey","black"))+
  theme_bw()

## Fig
ggplot(gg_qc,aes(x=seq.combined@active.ident,y=percent.mt,
                 alpha = 0.2, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_violin()+
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "grey","black"))+
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                             "darkorange3","coral3",
                             "darkred","grey","black"))+
  theme_bw()

ggplot(gg_qc,aes(x=nCount_RNA,
                 alpha = 0.05, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_density()+
  theme_bw()


#Setting cut-off values
dim(subset(seq.combined, 
           subset = 
             nFeature_RNA > 200 &
             nFeature_RNA < 6000 &
             nCount_RNA > 500 &
             nCount_RNA < 20000 &
             percent.mt < 10))

#Normalization
CP10K <- NormalizeData(seq.combined, normalization.method = "RC",
                       scale.factor = 10000)@assays$RNA
CP10K <- as.matrix(CP10K)
CP10K <- CP10K+0.1

seq.combined <- NormalizeData(seq.combined, normalization.method = "LogNormalize",
                        scale.factor = 10000)


#SC_RC9 <- FindVariableFeatures(seq.combined, selection.method = "vst", nfeatures = 2000)

###Seurat Tutorial###
seq.combined = FindVariableFeatures(seq.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(seq.combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seq.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the Data
all.genes = row.names(seq.combined)
seq.combined = ScaleData(seq.combined, features = all.genes)

seq.combined = SCTransform(seq.combined, method = "glmGamPoi",
                      vars.to.regress = c("nFeature_RNA","percent.mt","nonmt.libsize"),
                      return.only.var.genes = F)#"nFeature_RNA","percent.mt",

seq.combined <- RunPCA(seq.combined,
                 features = VariableFeatures(object = seq.combined))
eigValues <- (seq.combined@reductions$pca@stdev)**2
varExpl <- eigValues/sum(eigValues)

gg_pca <- data.frame(sample=row.names(seq.combined@reductions$pca@cell.embeddings),
                     PC1=seq.combined@reductions$pca@cell.embeddings[,1],
                     PC2=seq.combined@reductions$pca@cell.embeddings[,2],
                     PC3=seq.combined@reductions$pca@cell.embeddings[,3],
                     PC4=seq.combined@reductions$pca@cell.embeddings[,4],
                     PC5=seq.combined@reductions$pca@cell.embeddings[,5],
                     time=NA,
                     perc.mt=NA,
                     nonmt.libsize=seq.combined[["nonmt.libsize"]],
                     number.genes=seq.combined@meta.data$nFeature_RNA)

for(r in 1:dim(gg_pca)[1]){
  #gg_pca[r,7] <- as.character(tags[which(tags$X==gg_pca[r,1]),2])
  gg_pca[r,8] <- seq.combined[["percent.mt"]][which(row.names(seq.combined[["percent.mt"]])==gg_pca[r,1]),1]
}





## Fig 5.23 A
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))

#percent.mt
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=perc.mt, fill=perc.mt),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))
#nonmt.libsize ### DOES NOT WORK###
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=log10(nonmt.libsize), fill=log10(nonmt.libsize)),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))
#number.genes
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=log10(number.genes), fill=log10(number.genes)),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))

# PC2 vs PC3
ggplot(gg_pca)+
  geom_point(aes(x=PC2,y=PC3,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.2) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC2: ",round(varExpl[2]*100,1),"%"),
       y=paste0("PC3: ",round(varExpl[3]*100,1),"%"))
# PC3 vs PC4
ggplot(gg_pca)+
  geom_point(aes(x=PC3,y=PC4,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC3: ",round(varExpl[3]*100,1),"%"),
       y=paste0("PC4: ",round(varExpl[4]*100,1),"%"))

# PC4 vs PC5
ggplot(gg_pca)+
  geom_point(aes(x=PC4,y=PC5,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC4: ",round(varExpl[4]*100,1),"%"),
       y=paste0("PC5: ",round(varExpl[5]*100,1),"%"))

#Density vs PC2
ggplot(gg_pca,aes(x=PC2,
                  alpha = 0.05, color=seq.combined@active.ident, fill=seq.combined@active.ident)) +
  geom_density()+
  theme_bw()

#Centroids
mean(gg_pca[which(seq.combined@active.ident == "N2i_WT_1"),3])

samples = c("N2i_WT_1","N2i_WT_2","N12_WT_1","N12_WT_2","N24_WT_1","N24_WT_2")


centroid_frame <- data.frame("PC2"=rep(NA,6),
                             "PC3"=rep(NA,6),
                             row.names = samples)

for (r in 1:dim(centroid_frame)[1]){
  for (c in 1:dim(centroid_frame)[2]){
    centroid_frame[r,c] = mean(gg_pca[which(seq.combined@active.ident == row.names(centroid_frame)[r]),c+2])
  }
}

# Centroids PC2 vs PC3
ggplot(centroid_frame)+
  geom_point(aes(x=PC2,y=PC3,size=4, color=samples, fill=samples),
             size=2, alpha=1) +
  scale_color_manual(values=c("springgreen2","springgreen4",
                              "red2", "red4",
                              "darkgoldenrod1","darkgoldenrod3"
                              )) +
  theme_bw() +
  labs(x=paste0("PC2: ",round(varExpl[2]*100,1),"%"),
       y=paste0("PC3: ",round(varExpl[3]*100,1),"%"))


#Top Loadings
pca = RunPCA(seq.combined,
             features = VariableFeatures(object = seq.combined), 
             verbose = TRUE)
#How many cells of 24h are in the 0h time point -> 'not differentiated' -> 9.15%
count_of_value = sum(gg_pca$PC2[which(seq.combined@active.ident%in%c("N24_WT_1","N24_WT_2"))] <= mean(centroid_frame$PC2[1:2]))
percentage_of_value = (count_of_value/nrow(gg_pca[which(seq.combined@active.ident%in%c("N24_WT_1","N24_WT_2")),])) * 100
print(percentage_of_value)
#How many cells of 0h are in the 24h time point -> 'too fast differentiated' -> 0.54%
count_of_value = sum(gg_pca$PC2[which(seq.combined@active.ident%in%c("N2i_WT_1","N2i_WT_2"))] >= mean(centroid_frame$PC2[5:6]))
percentage_of_value = (count_of_value/nrow(gg_pca[which(seq.combined@active.ident%in%c("N2i_WT_1","N2i_WT_2")),])) * 100
print(percentage_of_value)

#UMAP Analysis
seq.combined = FindNeighbors(seq.combined, dims = 1:20, verbose = F)
seq.combined = FindClusters(seq.combined, resolution = 0.5, verbose = F)

seq.combined = RunUMAP(seq.combined, dims = 1:20, verbose = F)
gg_pca$UMAP_1 = seq.combined@reductions$umap@cell.embeddings[,1]
gg_pca$UMAP_2 = seq.combined@reductions$umap@cell.embeddings[,2]

ggplot(gg_pca)+
  geom_point(aes(x=UMAP_1,y=UMAP_2,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x="UMAP_1", y="UMAP_2")

#tSNE Analysis
seq.combined <- RunTSNE(seq.combined, dims = 1:20, verbose = FALSE)
gg_pca$tSNE_1 <- seq.combined@reductions$tsne@cell.embeddings[,1]
gg_pca$tSNE_2 <- seq.combined@reductions$tsne@cell.embeddings[,2]

ggplot(gg_pca)+
  geom_point(aes(x=tSNE_1,y=tSNE_2,size=4, color=seq.combined@active.ident, fill=seq.combined@active.ident),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x="tSNE_1", y="tSNE_2")

#Save seq.combined
saveRDS(seq.combined, file = "/cellfile/datapublic/fmosko/RWorkspace/seq.combined")

# retrieve a lookup between ensembl gene ids and mgi symbols from biomart
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
                 mart=ensembl)

sct_counts <- seq.combined@assays[["SCT"]]@counts
sct_counts <- sct_counts + 0.1

mgi_seq.combined <- rep(NA,nrow(sct_counts))
biotype_seq.combined <- rep(NA,nrow(sct_counts))
for(d in 1:dim(sct_counts)[1]){
  if(length(biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),2])==1){
    biotype_seq.combined[d] <- biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),2]
  }
  if(length(biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),2])>1){
    mart_types <- biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),2]
    if(length(unique(mart_types))==1){
      biotype_seq.combined[d] <- mart_types[1]
    }
    else{
      print(d)
    }
  }
  if(length(biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),2])==0){
    biotype_seq.combined[d] <- NA
  }
  if(length(biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),3])==TRUE){
    mgi_seq.combined[d] <- biomart[which(biomart$mgi_symbol==rownames(sct_counts)[d]),3]
  }
  else{
    mgi_seq.combined[d] <- NA
  }
}

##remove all genes with NA as mgi symbol
#seq.combined <- sct_counts[!is.na(mgi_seq.combined),]
biotype_seq.combined <- biotype_seq.combined[!is.na(mgi_seq.combined)]
mgi_seq.combined <- mgi_seq.combined[!is.na(mgi_seq.combined)]

##remove all genes with empty mgi symbol
sct_counts <- sct_counts[which(mgi_seq.combined!=""),]
biotype_seq.combined <- biotype_seq.combined[which(mgi_seq.combined!="")]
mgi_seq.combined <- mgi_seq.combined[which(mgi_seq.combined!="")]

##remove all genes with duplicate mgi symbols
sct_counts <- sct_counts[mgi_seq.combined%in%names(table(mgi_seq.combined)[which(table(mgi_seq.combined)==1)]),]
biotype_seq.combined <- biotype_seq.combined[mgi_seq.combined%in%names(table(mgi_seq.combined)[which(table(mgi_seq.combined)==1)])]
mgi_seq.combined <- mgi_seq.combined[mgi_seq.combined%in%names(table(mgi_seq.combined)[which(table(mgi_seq.combined)==1)])]

sct_counts <- sct_counts[which(biotype_seq.combined%in%
                       c("protein_coding","lincRNA",
                         "processed_transcript","antisense",
                         "3prime_overlapping_ncRNA",
                         "bidirectional_promoter_lncRNA",
                         "macro_lncRNA","miRNA",
                         "misc_RNA","lncRNA",
                         "scaRNA","scRNA",
                         "sense_intronic","sense_overlapping",
                         "snoRNA","snRNA","sRNA")),]
sct_counts <- sct_counts[which(row.names(sct_counts)%in%
                               row.names(sct_counts)),]

mgi_seq.combined <- mgi_seq.combined[which(biotype_seq.combined%in%
                               c("protein_coding","lincRNA",
                                 "processed_transcript","antisense",
                                 "3prime_overlapping_ncRNA",
                                 "bidirectional_promoter_lncRNA",
                                 "macro_lncRNA","miRNA",
                                 "misc_RNA","lncRNA",
                                 "scaRNA","scRNA",
                                 "sense_intronic","sense_overlapping",
                                 "snoRNA","snRNA","sRNA"))]
biotype_seq.combined <- biotype_seq.combined[which(biotype_seq.combined%in%
                                       c("protein_coding","lincRNA",
                                         "processed_transcript","antisense",
                                         "3prime_overlapping_ncRNA",
                                         "bidirectional_promoter_lncRNA",
                                         "macro_lncRNA","miRNA",
                                         "misc_RNA","lncRNA",
                                         "scaRNA","scRNA",
                                         "sense_intronic","sense_overlapping",
                                         "snoRNA","snRNA","sRNA"))]

## create gain loss frame
gain_loss_SC <- data.frame(TP_0 = rep(NA,nrow(sct_counts)),
                           TP_12 = rep(NA,nrow(sct_counts)),
                           TP_24 = rep(NA,nrow(sct_counts)),
                           change_sign_all = rep(NA,nrow(sct_counts)),
                           log2FC_avg = rep(NA,nrow(sct_counts)),
                           log2FC_12_0 = rep(NA,nrow(sct_counts)),
                           log2FC_24_0 = rep(NA,nrow(sct_counts)),
                           change_sign = rep(NA,nrow(sct_counts)),
                           row.names = row.names(sct_counts))

for(d in 1:nrow(gain_loss_SC)){
  if(d %% 1500 == 0){print(d)}
  ## T0
  gain_loss_SC[d,1] <- mean(as.numeric(sct_counts[d,which(seq.combined@active.ident%in%c("N2i_WT_1","N2i_WT_2"))]))
  ## T12
  gain_loss_SC[d,2] <- mean(as.numeric(sct_counts[d,which(seq.combined@active.ident%in%c("N12_WT_1","N12_WT_2"))]))
  ## T24
  gain_loss_SC[d,3] <- mean(as.numeric(sct_counts[d,which(seq.combined@active.ident%in%c("N24_WT_1","N24_WT_2"))]))
  
  ## is there a consistent gain in expression
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])%in%c(0,1) &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])%in%c(0,1) #&
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])%in%c(0,1)
  ){
    gain_loss_SC[d,4] <- 1
  }
  ## is there a consistent loss in expression
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])%in%c(0,-1) &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])%in%c(0,-1) #&
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])%in%c(0,-1)
  ){
    gain_loss_SC[d,4] <- -1
  }
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])==0 &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])==0 
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])==0
  ){
    gain_loss_SC[d,4] <- 0
  }
  ## log2 FC all others vs T0
  gain_loss_SC[d,5] <- mean(c(log2((gain_loss_SC[d,2])/(gain_loss_SC[d,1])),
                              log2((gain_loss_SC[d,3])/(gain_loss_SC[d,1]))))
  ## log2 FC all T12 vs T0
  gain_loss_SC[d,6] <- log2((gain_loss_SC[d,2])/(gain_loss_SC[d,1]))
  ## log2 FC all T24 vs T0
  gain_loss_SC[d,7] <- log2((gain_loss_SC[d,3])/(gain_loss_SC[d,1]))
  
}

gain_loss_ggplot <- gain_loss_SC[,4:7]
gain_loss_ggplot <- gain_loss_ggplot[rowSums(is.na(gain_loss_ggplot))!=ncol(gain_loss_ggplot),]
## plot log2 FC T24 vs T0 against average log2 FC vs T0
ggplot(gain_loss_ggplot, aes(x=log2FC_avg, y=log2FC_24_0,
                             color=as.factor(change_sign_all),
                             fill=as.factor(change_sign_all))) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              col="red") +
  labs(x="average log2 FC vs 0h",
       y="log2 FC 24h vs 0h")
## plot log2 FC T48 vs T0 against log2 FC T24 vs T0
ggplot(gain_loss_ggplot, aes(x=log2FC_12_0, y=log2FC_24_0,
                             color=as.factor(change_sign_all),
                             fill=as.factor(change_sign_all))) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              col="red") +
  labs(x="log2 FC 12h vs 0h",
       y="log2 FC 24h vs 0h")

saveRDS(seq.combined, "/cellfile/datapublic/fmosko/RWorkspace/seq.combined.rds")
saveRDS(sct_counts, "/cellfile/datapublic/fmosko/RWorkspace/sct_counts.rds")
saveRDS(gain_loss_SC, "/cellfile/datapublic/fmosko/RWorkspace/gain_loss_SC.rds")

saveRDS(data.frame(ensembl=row.names(seq.combined),
                  mgi=mgi_seq.combined),"RDS/log2_seq.combined_lookup.rds")

### start script 08 ###

seq.combined = readRDS("/cellfile/datapublic/fmosko/RWorkspace/seq.combined.rds")
sct_counts = readRDS("/cellfile/datapublic/fmosko/RWorkspace/sct_counts.rds")
gain_loss_SC = readRDS("/cellfile/datapublic/fmosko/RWorkspace/gain_loss_SC.rds")

colnames(sct_counts) <- paste(seq.combined@active.ident,
                           1:length(seq.combined@active.ident),sep="_")
sc_data_0to24 <- sct_counts

T0_cells <- which(seq.combined@active.ident%in%c("N2i_WT_1","N2i_WT_2"))
T12_cells <- which(seq.combined@active.ident%in%c("N12_WT_1","N12_WT_2"))
T24_cells <- which(seq.combined@active.ident%in%c("N24_WT_1","N24_WT_2"))

time_points_0to24 <- sapply(colnames(sc_data_0to24),
                            function(x) strsplit(x,split="_")[[1]][1])

### subsamplen ###
sample_cells <- 1:length(colnames(sc_data_0to24))
names(sample_cells) <- colnames(sc_data_0to24)
samples = c("N2i_WT_1", "N2i_WT_2", "N12_WT_1", "N12_WT_2", "N24_WT_1", "N24_WT_2")

sample_numbers <- data.frame("Numbers"=rep(NA,6),
                              row.names = samples)

for (r in 1:dim(sample_numbers)[1]){
  sample_numbers[r, "Numbers"] = sum(seq.combined@active.ident == row.names(sample_numbers)[r])
  
}

set.seed(123)
for (i in samples) {
  sample_cells[grep(i,names(sample_cells))] = sample(1:sample_numbers[i,1], sample_numbers[i,1], replace = FALSE)
}

sample_cells[grep(5014,sample_cells)]
#sample_cells_b1 = sample_cells[which(sample_cells <= 1000)]
#sample_cells_b2 = sample_cells[which(sample_cells > 1000 & sample_cells <= 2000)]

### subset with 1000 cells per sample
sc_exp_obj <- SingleCellExperiment(list(count=sc_data_0to24[,which(sample_cells > 3000 & sample_cells <= 6000)]))
logcounts(sc_exp_obj) <- log2(sc_data_0to24[,which(sample_cells > 3000 & sample_cells <= 6000)])

#saveRDS(sc_exp_obj, "/cellfile/datapublic/fmosko/RWorkspace/sc_exp_obj.rds")
#saveRDS(time_points_0to24, "/cellfile/datapublic/fmosko/RWorkspace/time_points_0to24.rds")

Sys.time()
pseudo_obj_b2 <- psupertime(x = sc_exp_obj,
                         y = factor(time_points_0to24[which(sample_cells > 3000 & sample_cells <= 6000)],
                                    levels = c("N2i", "N12", "N24")),
                         sel_genes = "all")
Sys.time()

#saveRDS(pseudo_obj, "/cellfile/datapublic/fmosko/RWorkspace/pseudo_obj.rds")

plot_train_results(pseudo_obj_b2)
plot_labels_over_psupertime(pseudo_obj_b2, label_name='Time', palette = "YlOrRd")
plot_identified_gene_coefficients(pseudo_obj_b2)
plot_identified_genes_over_psupertime(pseudo_obj_b1, label_name='Time', palette = "YlOrRd")

saveRDS(pseudo_obj_b1, "/cellfile/datapublic/fmosko/RWorkspace/pseudo_obj_b1.rds")
saveRDS(pseudo_obj_b2, "/cellfile/datapublic/fmosko/RWorkspace/pseudo_obj_b2.rds")
saveRDS(names_lookup, "/cellfile/datapublic/fmosko/RWorkspace/names_lookup.rds")
saveRDS(sample_numbers, "/cellfile/datapublic/fmosko/RWorkspace/sample_numbers.rds")
saveRDS(sc_exp_obj, "/cellfile/datapublic/fmosko/RWorkspace/sc_exp_obj.rds")
saveRDS(sc_data_0to24, "/cellfile/datapublic/fmosko/RWorkspace/sc_data_0to24.rds")
saveRDS(sc_data_0to24_ordered_b1, "/cellfile/datapublic/fmosko/RWorkspace/sc_data_0to24_ordered_b1.rds")
saveRDS(sc_data_0to24_ordered_b2, "/cellfile/datapublic/fmosko/RWorkspace/sc_data_0to24_ordered_b2.rds")
saveRDS(genes_summary, "/cellfile/datapublic/fmosko/RWorkspace/genes_summary.rds")

#batch1 vs batch2! Where are the top20 genes of batch1 in batch2?
which(pseudo_obj_b1[["beta_dt"]][["symbol"]] %in% pseudo_obj_b2[["beta_dt"]][["symbol"]][1:20])
max(which(pseudo_obj_b1[["beta_dt"]][["symbol"]] %in% pseudo_obj_b2[["beta_dt"]][["symbol"]][1:20]))
#And vice-versa
which(pseudo_obj_b2[["beta_dt"]][["symbol"]] %in% pseudo_obj_b1[["beta_dt"]][["symbol"]][1:20])
max(which(pseudo_obj_b2[["beta_dt"]][["symbol"]] %in% pseudo_obj_b1[["beta_dt"]][["symbol"]][1:20]))


pseudo_order_b2 <- order(pseudo_obj_b2$proj_dt$psuper, decreasing = F)
sc_data_0to24_ordered_b2 <- sc_data_0to24[,pseudo_order_b2]



## biomart lookup between ensembl and mgi symbols
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")#, host = "uswest.ensembl.org")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart),]
### cell cycle prediction ####
## retrieve a lookup table for ensembl and mgi gene names
biomart_cc <- biomart[which(!is.na(biomart$entrezgene_id)),]
biomart_cc <- biomart_cc[complete.cases(biomart_cc),]

names_lookup <- data.frame(mgi=sapply(rownames(sc_data_0to24),function(x) strsplit(x,"_")[[1]][1]),
                           ensembl=NA)
for(d in 1:dim(names_lookup)[1]){
  if(length(biomart_cc[which(biomart_cc$mgi_symbol==names_lookup[d,1]),1])==1){
    names_lookup[d,2] <- biomart_cc[which(biomart_cc$mgi_symbol==names_lookup[d,1]),1][1]
  }
}

## sort out esembl genes that were not unique
unique_ensembl <- names_lookup[which(!names_lookup$ensembl%in%
                                       unique(names_lookup[duplicated(names_lookup[,2]),2])),2]

## which genes are unique in the ensembl lookup and can be used in the cellcycle
## calculation
#sc_data_0to24_cellcycle <- sc_data_0to24[which(row.names(sc_data_0to24)%in%unique_ensembl),]

#mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
#cell_cycle <- cyclone(as.matrix(sc_data_0to24_cellcycle),mm.pairs)
#names(cell_cycle$phases) <- colnames(sc_data_0to24)


### identify genes that expressed at certain level in most cells ####
A_relation_b_hm_consist <- readRDS("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/RDS/A_relation_b_hm_consist_perc_genes.rds")
## set rownames and colnames to ensembl for easy handling
A_relation_b_hm_consist <- A_relation_b_hm_consist[row.names(A_relation_b_hm_consist)%in%na.omit(names_lookup)[,1],
                                                   colnames(A_relation_b_hm_consist)%in%na.omit(names_lookup)[,1]]
for(d in 1:nrow(A_relation_b_hm_consist)){
  row.names(A_relation_b_hm_consist)[d] <- names_lookup[which(names_lookup$mgi==row.names(A_relation_b_hm_consist)[d]),2]
  colnames(A_relation_b_hm_consist)[d] <- row.names(A_relation_b_hm_consist)[d]
}

hist(unlist(as.matrix(log2(sct_counts)),use.names = F))

## summarize how big the percentage of cells is that make a certain TPM
## cutoff either for all cells or for certain time points.
summary_thres <- data.frame(gene=row.names(sct_counts),
                            ## columns summarizing the percentage of cells making the cut off over all
                            ## cells
                            bigger_1=apply(sct_counts,1,function(x) (sum(x>1,na.rm=T)/
                                                                    length(x))),
                            bigger_2=apply(sct_counts,1,function(x) (sum(x>2,na.rm=T)/
                                                                    length(x))),
                            bigger_3=apply(sct_counts,1,function(x) (sum(x>3,na.rm=T)/
                                                                    length(x))),
                            bigger_4=apply(sct_counts,1,function(x) (sum(x>4,na.rm=T)/
                                                                    length(x))),
                            bigger_5=apply(sct_counts,1,function(x) (sum(x>5,na.rm=T)/
                                                                    length(x))),
                            ##  the percentage of cells having NAs over all cells
                            nas=apply(sct_counts,1,function(x) (sum(x<0.1)/length(x))),#(sum(is.na(x))/length(x))),
                            ## columns summarizing the percentage of cells making the cut off over all
                            ## cells per time point
                            bigger_N2i_WT_1=apply(sct_counts[,which(seq.combined@active.ident=="N2i_WT_1")]
                                              ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_N2i_WT_2=apply(sct_counts[,which(seq.combined@active.ident=="N2i_WT_2")]
                                              ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_N12_WT_1=apply(sct_counts[,which(seq.combined@active.ident=="N12_WT_1")]
                                              ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_N12_WT_2=apply(sct_counts[,which(seq.combined@active.ident=="N12_WT_2")]
                                              ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_N24_WT_1=apply(sct_counts[,which(seq.combined@active.ident=="N24_WT_1")]
                                               ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_N24_WT_2=apply(sct_counts[,which(seq.combined@active.ident=="N24_WT_2")]
                                               ,1,function(x) (sum(x>2,na.rm=T)/length(x))))
## prepare plot to compare cut offs
gg_thresh <- data.frame(percentage=c(summary_thres$bigger_1,
                                     summary_thres$bigger_2,
                                     summary_thres$bigger_3),
                        group=c(rep("count > 1",length(summary_thres$bigger_1)),
                                rep("count > 2",length(summary_thres$bigger_2)),
                                rep("count > 3",length(summary_thres$bigger_3))))

ggplot(gg_thresh, aes(x=percentage,color=group,fill=group)) +
  geom_density(alpha=0.2)+
  ggtitle("density plot of genes passing different expression thresholds") +
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(0,0.05)

## read in time course data and define which genes change in the TC
gpr_list_shrunkenFCS <- readRDS("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/RDS/gpr_list_shrunkenFCS.rds")
log_norm_heatmap_gpr <- gpr_list_shrunkenFCS[[1]][,2:18]
log_norm_heatmap_gpr <- log_norm_heatmap_gpr-log_norm_heatmap_gpr[,1]
log_norm_heatmap_gpr_changed <- log_norm_heatmap_gpr[which(rowMax(log_norm_heatmap_gpr)-
                                                             rowMin(log_norm_heatmap_gpr)>=0.5),]
changed_genes_ensembl <- row.names(log_norm_heatmap_gpr_changed)
changed_genes_mgi <- biomart[which(biomart$ensembl_gene_id%in%changed_genes_ensembl),3]



## further summary what genes have a certain expression in at least
## 70 or 80 percent of the cells
genes_summary <- data.frame(row.names=row.names(sct_counts),
                            gene=names_lookup$mgi,
                            ## columns summarizing whether a gene has at least 50 or 60
                            ## percent of cells above a TPM cutoff
                            TPM_bigger_1_20_perc=summary_thres$bigger_1>=0.2,
                            TPM_bigger_1_30_perc=summary_thres$bigger_1>=0.3,
                            ## columns summarizing whether a gene has less than 10
                            ## percent of NAs per gene, mean, median and variance TPM per gene
                            Na_perc_lower_10=summary_thres$nas<0.1,
                            mean_TPM=apply(sct_counts,1,function(x) mean(x,na.rm=T)),
                            median_TPM=apply(sct_counts,1,function(x) median(x,na.rm=T)),
                            var_TPM=apply(sct_counts,1,function(x) var(x,na.rm=T)),
                            ## columns summarizing whether a gene has at least 70 or 80
                            ## percent of cells above a TPM cutoff per time point
                            TPM_bigger_1_20_perc_N2i_WT_1=summary_thres$bigger_N2i_WT_1>=0.2,
                            TPM_bigger_1_30_perc_N2i_WT_1=summary_thres$bigger_N2i_WT_1>=0.3,
                            TPM_bigger_1_20_perc_N2i_WT_2=summary_thres$bigger_N2i_WT_2>=0.2,
                            TPM_bigger_1_30_perc_N2i_WT_2=summary_thres$bigger_N2i_WT_2>=0.3,
                            TPM_bigger_1_20_perc_N12_WT_1=summary_thres$bigger_N12_WT_1>=0.2,
                            TPM_bigger_1_30_perc_N12_WT_1=summary_thres$bigger_N12_WT_1>=0.3,
                            TPM_bigger_1_20_perc_N12_WT_2=summary_thres$bigger_N12_WT_2>=0.2,
                            TPM_bigger_1_30_perc_N12_WT_2=summary_thres$bigger_N12_WT_2>=0.3,
                            TPM_bigger_1_20_perc_N24_WT_1=summary_thres$bigger_N24_WT_1>=0.2,
                            TPM_bigger_1_30_perc_N24_WT_1=summary_thres$bigger_N24_WT_1>=0.3,
                            TPM_bigger_1_20_perc_N24_WT_2=summary_thres$bigger_N24_WT_2>=0.2,
                            TPM_bigger_1_30_perc_N24_WT_2=summary_thres$bigger_N24_WT_2>=0.3)

## for every gene check where highest and lowest expression is in TC (window)
## then give each gene a sign based on that + high expression late, - low expression late
gpr_TPMs <- readRDS("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/RDS/tpms_TC_gpr.rds")
gpr_TPMs_short <- data.frame(RC9_2i=rowMeans(gpr_TPMs[,1:2]),
                             RC9_2h=gpr_TPMs[,5],
                             RC9_4h=gpr_TPMs[,6],
                             RC9_6h=gpr_TPMs[,7],
                             RC9_8h=gpr_TPMs[,8],
                             RC9_10h=gpr_TPMs[,9],
                             RC9_12h=gpr_TPMs[,10],
                             RC9_14h=gpr_TPMs[,11],
                             RC9_16h=gpr_TPMs[,12],
                             RC9_18h=gpr_TPMs[,13],
                             RC9_20h=gpr_TPMs[,14],
                             RC9_22h=gpr_TPMs[,15],
                             RC9_24h=gpr_TPMs[,16],
                             RC9_26h=gpr_TPMs[,17],
                             RC9_28h=gpr_TPMs[,18],
                             RC9_30h=gpr_TPMs[,19],
                             RC9_32h=rowMeans(gpr_TPMs[,20:21]),
                             mgi=gpr_TPMs[,23],
                             min_TP=NA,
                             max_TP=NA,
                             type=NA)

## use windows (3 time points together) to be less vulnerable to a single outlier
for(d in 1:nrow(gpr_TPMs_short)){
  windows <- c(rowMeans(gpr_TPMs_short[d,1:3]),
               rowMeans(gpr_TPMs_short[d,2:4]),
               rowMeans(gpr_TPMs_short[d,3:5]),
               rowMeans(gpr_TPMs_short[d,4:6]),
               rowMeans(gpr_TPMs_short[d,5:7]),
               rowMeans(gpr_TPMs_short[d,6:8]),
               rowMeans(gpr_TPMs_short[d,7:9]),
               rowMeans(gpr_TPMs_short[d,8:10]),
               rowMeans(gpr_TPMs_short[d,9:11]),
               rowMeans(gpr_TPMs_short[d,10:12]),
               rowMeans(gpr_TPMs_short[d,11:13]),
               rowMeans(gpr_TPMs_short[d,12:14]),
               rowMeans(gpr_TPMs_short[d,13:15]),
               rowMeans(gpr_TPMs_short[d,14:16]),
               rowMeans(gpr_TPMs_short[d,15:17]))
  ## if all windows are min and max there is no change and we
  ## fill NA as min and max
  if(length(which(windows==min(windows)))==15 &&
     length(which(windows==max(windows)))==15){
    gpr_TPMs_short[d,19] <- NA
    gpr_TPMs_short[d,20] <- NA
    ## fill in min and max of gene window per gene
  } else{
    gpr_TPMs_short[d,19] <- which(windows==min(windows))[1]
    gpr_TPMs_short[d,20] <- which(windows==max(windows))[1]
    ## if window min is located before window max sign is positive
    if(gpr_TPMs_short[d,19] < gpr_TPMs_short[d,20]){
      gpr_TPMs_short[d,21] <- 1
      ## if window min is located after window max sign is negative
    } else if(gpr_TPMs_short[d,19] > gpr_TPMs_short[d,20]){
      gpr_TPMs_short[d,21] <- -1
      ## else no sign
    } else{
      gpr_TPMs_short[d,21] <- 0
    }
    
  }
}

## remove the directions that are inconsistent with SC data
## here we look at the changes between E5.5 and E3.5 as E4.5 has
## a very small number of cells
gain_loss_ggplot$change_24_vs_0 <- NA
for(d in 1:nrow(gain_loss_ggplot)){
  ## log2FCs E4.5 vs E3.5 and E5.5 vs E3.5 should not be na
  if(!is.na(gain_loss_ggplot[d,5])){
    ## if log2FC E5.5 vs E3.5 is positive sign is positive
    if(sign(gain_loss_ggplot[d,5])==1){
      gain_loss_ggplot[d,7] <- 1
    }
    ## if log2FC E5.5 vs E3.5 is negative sign is negative
    if(sign(gain_loss_ggplot[d,5])==-1){
      gain_loss_ggplot[d,7] <- -1
    }
  }
}

saveRDS(gpr_TPMs_short, "/cellfile/datapublic/fmosko/RWorkspace/genes_summary.rds")
