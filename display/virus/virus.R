#!/usr/bin/env Rscript
suppressMessages({library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)})
#readRDS
rds=readRDS("/SGRNJ02/RandD4/Integratedanalysis/virus_panel/20210107/amulti/virus.diff_PRO.rds")
df=read_tsv("/SGRNJ02/RandD4/virus_panel/20210102/S20201218_EBV__EBV_KZ/05.analysis_capture_virus/S20201218_EBV__EBV_KZ_tsne.tsv")
#do_go
do_go <- function(genes){
    converted <- select(org.Hs.eg.db, 
       keys = genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

    gene_id <- na.omit(converted$ENTREZID)

    ego <- enrichGO(gene          = gene_id,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
    return (ego)
}
####
df$shi<-c('S20201218_EBV__EBV_')
df$barcode<-str_c(df$shi,df$barcode)
df<-df[,-9]
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
####
meta = rds@meta.data
meta$vir = 'Not_detected'
df.vir = df[grepl('EBV', df$tag),]
df.vir=df.vir[df.vir$UMI>5,]
vir.barcode = df.vir$barcode
meta[vir.barcode,'vir'] = 'EBV'
rds@meta.data=meta
new.cluster.ids <- c('Tcells','Tcells','Monocytes','Bcells','Tcells','Endothelial cells','Plasma cells')
names(new.cluster.ids) <- levels(rds)
rds <- RenameIdents(rds, new.cluster.ids)
rds <- StashIdent(object = rds, save.name = "ClusterNames")
P<-TSNEPlot(object = rds,label = TRUE)
print(P)
ggplot2::ggsave(filename = "S20201218_EBV_EBV_lab_tsne.pdf")
Idents(rds) <- "vir"
P<-TSNEPlot(object = rds)
print(P)
ggplot2::ggsave(filename = "S20201218_EBV_EBV_vir_tsne.pdf")
#######
markergene=c('RAG1', 'RAG2', 'DCLRE1C', 'LIG4','IL2RG', 'JAK3', 'IL7RA','STIM1', 'RASGRP1' , 'CTPS1' , 'MST1', 'GATA2', 'DOCK8', 'WAS','CORO1A', 'PIK3CD‐GoF', 'PIK3R1', 'NFKB1','ITK', 'MAGT1','CD27', 'CD70', 'TNFRSF9')
Idents(rds) <- "vir"
cluster.markers <- FindMarkers(object=rds,ident.1="EBV",ident.2='Not_detected',logfc.threshold = 0)
write.table(cluster.markers,file="EBV_diffgenes.xls",,sep='\t',quote=F,row.names=T)
#IFITM3<-FindMarkers(object=rds,ident.1="EBV",ident.2='Not_detected',features = markergene,logfc.threshold = 0)
#write.table(IFITM3,file="EBV_IFITM3diffgenes.xls",,sep='\t',quote=F,row.names=T)
pdf(file = "markergene.pdf",width =12,height = 15)
VlnPlot(object = rds, features =markergene)
dev.off()
rds<-FindVariableFeatures(rds)
genes.use <- head(HVFInfo(object = rds))
all_data.markers <- FindAllMarkers(object = rds, genes.use = genes.use)
write.csv(all_data.markers,"all_markers.csv")
top_10 <- all_data.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top_10,"top_10_markers.csv")
for (c_cluster in unique(all_data.markers$cluster)){
    c_markers <- filter(all_data.markers,cluster==c_cluster)
    file_name <- paste0(c_cluster,"_markers.tsv")
    write_tsv(c_markers,file_name)
}
for (cluster in unique(all_data.markers$cluster)){
  markers <- as.vector(all_data.markers[all_data.markers$cluster==cluster & all_data.markers$avg_logFC >0,]$gene)
  tryCatch({
    ego <- do_go(markers)

    result_file <- paste(cluster,"_go.csv",sep="")
    write.csv(ego@result,result_file)
    
    pdf(file = paste(cluster,"_go_dot.pdf",sep=""),width =15,height = 12)
    print (dotplot(ego,showCategory=20))
    dev.off()


    pdf(file = paste(cluster,"_go_bar.png",sep=""),width =15,height = 12)
    print (barplot(ego,showCategory=20))
    dev.off()
  
    },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})
}
Idents(rds) <- "ClusterNames"
Tcells<-subset(rds, idents = c("Tcells"))
Idents(Tcells) <- "vir"
cluster.markers <- FindMarkers(object=Tcells,ident.1="EBV",ident.2='Not_detected',logfc.threshold = 0)
write.table(cluster.markers,file="Tcells_diffgenes.xls",,sep='\t',quote=F,row.names=T)
pdf(file = "Tcells_markergene.pdf",width =12,height = 15)
VlnPlot(object = Tcells, features =markergene )
dev.off()
Tcells<-FindVariableFeatures(Tcells)
genes.use <- head(HVFInfo(object = Tcells))
Tcells_data.markers <- FindAllMarkers(object = Tcells, genes.use = genes.use)
write.csv(Tcells_data.markers,"Tcells_markers.csv")
top_10 <- Tcells_data.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top_10,"Tcells_top_10_markers.csv")

for (c_cluster in unique(Tcells_data.markers$cluster)){
    c_markers <- filter(Tcells_data.markers,cluster==c_cluster)
    file_name <- paste0("Tcells",c_cluster,"_markers.tsv")
    write_tsv(c_markers,file_name)
}

for (cluster in unique(Tcells_data.markers$cluster)){
  markers <- as.vector(Tcells_data.markers[Tcells_data.markers$cluster==cluster & Tcells_data.markers$avg_log2FC >0,]$gene)
  tryCatch({
    ego <- do_go(markers)
    result_file <- paste("Tcells",cluster,"_go.csv",sep="")
    write.csv(ego@result,result_file)
    pdf(file = paste("Tcells",cluster,"_go_dot.pdf",sep=""),width =15,height = 12)
    print (dotplot(ego,showCategory=20))
    dev.off()
    pdf(file = paste("Tcells",cluster,"_go_bar.pdf",sep=""),width =15,height = 12)
    print (barplot(ego,showCategory=20))
    dev.off()
    },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})
}
########
rds.vir<-rds@meta.data[grepl('EBV', rds@meta.data$vir),]
freq_table <- prop.table(x=table(rds.vir$ClusterNames,rds.vir[,"orig.ident"]),margin=2)
write.table(freq_table,file='freq_table.xls',sep='\t',quote=F,row.names=T,col.names = NA)
barplot(height=freq_table,width = 20,xlim=c(1,60),col =clustcol,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
pdf(file = "percentage.pdf",width =12,height = 15)
barplot(height=freq_table,width = 20,xlim=c(1,60),col =clustcol,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
dev.off()
