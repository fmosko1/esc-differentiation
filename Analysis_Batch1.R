gain_loss_SC = readRDS("/cellfile/datapublic/fmosko/RWorkspace/gain_loss_SC.rds")
gpr_TPMs_short = readRDS("/cellfile/datapublic/fmosko/RWorkspace/gpr_TMPs_short.rds")

gain_loss_ggplot <- gain_loss_SC[,4:7]
gain_loss_ggplot <- gain_loss_ggplot[rowSums(is.na(gain_loss_ggplot))!=ncol(gain_loss_ggplot),]

genes_one_direct <- row.names(gain_loss_ggplot[which(gain_loss_ggplot$change_sign_all%in%c(-1,1)),])


## prepare a plot to summarize log2FCs in SC data
gg_SC_TC <- gain_loss_ggplot[,c(3,4)]
gg_SC_TC$TC_sign <- NA
for(d in 1:nrow(gg_SC_TC)){
  if(row.names(gg_SC_TC)[d]%in%gpr_TPMs_short$mgi){
    if(length(gpr_TPMs_short[which(gpr_TPMs_short$mgi==row.names(gg_SC_TC)[d]),21])==1){
      gg_SC_TC[d,3] <- gpr_TPMs_short[which(gpr_TPMs_short$mgi==row.names(gg_SC_TC)[d]),21]
    }
  }
}

#Comparison of scSeq and bulk time course data
ggplot(gg_SC_TC, aes(x=log2FC_12_0,y=log2FC_24_0,
                     col=as.factor(TC_sign),fill=as.factor(TC_sign)))+
  geom_point()+
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  theme_bw() +
  geom_hline(yintercept = 0,
              col="black") +
  geom_vline(xintercept = 0,
             col="black") +
  labs(x="log2FC 12 vs 0",
       y="log2FC 24 vs 0") +
  ggtitle("consisctency change in SC data") +
  theme(plot.title=element_text(hjust=0.5))

# are directions of change between time course and SC data the same
gpr_TPMs_short_consist <- gpr_TPMs_short
for(d in 1:nrow(gpr_TPMs_short_consist)){
  if(gpr_TPMs_short_consist[d,21]%in%c(1,-1)){
    if(row.names(gpr_TPMs_short_consist)[d] %in% row.names(gain_loss_ggplot)){
      if(sign(gain_loss_ggplot[row.names(gpr_TPMs_short_consist)[d],5])!=gpr_TPMs_short_consist[d,21]){
        gpr_TPMs_short_consist[d,21] <- 0
      }
    } else {
      gpr_TPMs_short_consist[d,21] <- 0
    }
  }
}

#changed_genes_signs <- changed_genes_ensembl[which(changed_genes_ensembl %in%
#                                                     row.names(gpr_TPMs_short_consist[which(gpr_TPMs_short_consist[,21]%in%c(-1,1)),]))]
genes_one_direct_consist <- genes_one_direct[which(genes_one_direct %in%
                                                     row.names(gpr_TPMs_short_consist[which(gpr_TPMs_short_consist[,21]%in%c(-1,1)),]))]

### get thresholds and min max from the rocs ####
gpr_TPMs_short_consist$thresh <- NA
gpr_TPMs_short_consist$roc_FPR <- NA
gpr_TPMs_short_consist$roc_TPR <- NA
gpr_TPMs_short_consist$max_roc_FPR <- NA
gpr_TPMs_short_consist$min_roc_TPR <- NA

T0_cells <- which(time_points_0to24=="0h")
T24_cells <- which(time_points_0to24=="24h")


for(d in 1:nrow(gpr_TPMs_short_consist)){
  if(d%%1000 == 0){
    print(d)
  }
  if(gpr_TPMs_short_consist[d,21]%in%c(1,-1)){ #&gpr_TPMs_short_consist[d,18]%in%genes_pass
    if(sum(is.na(sc_data_0to24[row.names(gpr_TPMs_short_consist)[d],
                               c(T0_cells,T24_cells)]))==0){
      holder <- ROC_gg(gene=as.character(row.names(gpr_TPMs_short_consist)[d]),
                       SC_frame=sc_data_0to24, #log2_CP10k_0to24
                       direction=gpr_TPMs_short_consist[d,21],
                       nbreaks = 30, mode="thresh",
                       early = T0_cells, late = T24_cells)
      gpr_TPMs_short_consist[d,22] <- holder[1]
      gpr_TPMs_short_consist[d,23] <- holder[2]
      gpr_TPMs_short_consist[d,24] <- holder[3]
      
      holder <- ROC_gg(gene=as.character(row.names(gpr_TPMs_short_consist)[d]),
                       SC_frame=sc_data_0to24,
                       direction=gpr_TPMs_short_consist[d,21],
                       nbreaks = 30, mode="thresh",
                       NA_cut = 0.1,
                       early = T0_cells, late = T24_cells, min_TPR=T)
      gpr_TPMs_short_consist[d,25] <- holder[2]
      gpr_TPMs_short_consist[d,26] <- holder[3]
    }
    
  }
}
print("hello world")