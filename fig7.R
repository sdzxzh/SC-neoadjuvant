library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

#fig7a
#refer to fig4b

#fig7b
#refer to fig4c

#fig7c
library("Scillus")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
fib.cell <- rownames(subset(data.combined@meta.data, origin == "normal"))
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined1=sub_CT
data.combined1$fincell=Idents(data.combined1)
sce=data.combined1
sce$seurat_clusters=1
Idents(sce)<-"seurat_clusters"
sce$neo_group=sce$group
df13 <- find_diff_genes(dataset = sce, 
                      clusters = as.character(1),
                      comparison = c("group", "A", "D"),
                      logfc.threshold = 0,   
                      min.cells.group = 1)   
fib.cell <- rownames(subset(data.combined@meta.data, origin == "tumor"))
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined1=sub_CT
data.combined1$fincell=Idents(data.combined1)

library("Scillus")
sce=data.combined1
sce$seurat_clusters=1
Idents(sce)<-"seurat_clusters"
sce$neo_group=sce$group
df19 <- find_diff_genes(dataset = sce, 
                        clusters = as.character(1),
                        comparison = c("group", "A", "D"),
                        logfc.threshold = 0,   
                        min.cells.group = 1)   

rownames(df13)=df13$feature
rownames(df19)=df19$feature
gene=intersect(df13$feature,df19$feature)
df=cbind(df13[gene,],df19[gene,])
rownames(df)=df$feature
genes.to.label3 = intersect(genes.to.label1,genes.to.label2)
genes.to.label_fin=c(genes.to.label1[1:30],genes.to.label2[1:30],genes.to.label3)
genes.to.label_fin=genes.to.label_fin[-which(genes.to.label_fin%in%vgene)]
genes_to_showname=genes.to.label_fin
data=data.frame(df[,4],df[11])
rownames(data)=rownames(df)
colnames(data)=c("x","y")
data$gene_label <- ""
data[genes_to_showname, ]$gene_label <- genes_to_showname
xmin <- -0.25    
xmax <- 0.25     
ymin <- -0.25 
ymax <- 0.25  
x_axis_min <- -3  
x_axis_max <- 3   
y_axis_min <- -3  
y_axis_max <- 3   

xlab <- "log2(ratio of mRNA)"   
ylab <- "log2(ratio of protein)"   

x_tick_pos <- seq(x_axis_min, x_axis_max, 0.5)
y_tick_pos <- seq(y_axis_min, y_axis_max, 0.5)
point_color <- c("#d002ea", "#ff7c00", "#FB7D80", "grey70", "grey70", "#ff7c00", "#AF1E1F", "grey70", "grey70") 
names(point_color) <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9")
data <- data %>% mutate(quadrant = case_when(
  x > xmax & y > ymax ~ "Q1",
  x > xmax & y <= ymax & y > ymin ~ "Q2",
  x > xmax & y <= ymin ~ "Q3",
  x <= xmax & x > xmin & y <= ymin ~ "Q4",
  x <= xmax & x > xmin & y > ymin & y <= ymax ~ "Q5",
  x <= xmax & x > xmin & y > ymax ~ "Q6",
  x <= xmin & y > ymax ~ "Q7",
  x <= xmin & y <= ymax & y > ymin ~ "Q8",
  x <= xmin & y <= ymin ~ "Q9"))
head(data)

annotate_x <- c(rep(x_axis_max, 3), rep((x_axis_min+x_axis_max)/2, 3), rep(x_axis_min, 3))
annotate_y <- rep(c(y_axis_max, (y_axis_min+y_axis_max)/2, y_axis_min), 3)
annotate_text_color <- c(rep("black", 4), "white", rep("black", 4)) 

p <- ggplot(data) + 
  geom_point(aes(x=x, y=y, color=quadrant), size=1.5,alpha=0.6) + 
  coord_cartesian(xlim=c(x_axis_min, x_axis_max), 
                  ylim=c(y_axis_min, y_axis_max)) + 
  geom_vline(xintercept=c(xmin, xmax), size=0.3, linetype="dashed") + 
  geom_hline(yintercept=c(ymin, ymax), size=0.3, linetype="dashed") +
  annotate("text", x=annotate_x, y=annotate_y, label=c(1,2,3,6,5,4,7,8,9), color=annotate_text_color) + 
  labs(x=xlab, y=ylab) + 
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank()) + 
  scale_colour_manual(values=point_color) + 
  scale_x_continuous(breaks=x_tick_pos) +
  scale_y_continuous(breaks=y_tick_pos)
if (length(genes_to_showname) > 0) {
  p <- p + 
    geom_text_repel(aes(x=x, y=y, label=gene_label, color = quadrant), 
                    size=3, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"))
}
hist_x <- ggplot(data, aes(x=x)) + 
  geom_histogram(bins=100, 
                 color="black", fill = "#FB7D80",alpha=1) +
  coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  xlab(NULL)
hist_y <- ggplot(data, aes(x=y)) + 
  geom_histogram(bins=100, 
                 color="black", fill = "#AF1E1F",alpha=1) +
  coord_flip(xlim=c(y_axis_min, y_axis_max)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() ,axis.title.y=element_blank()) +
  xlab(NULL)
empty <- ggplot() + theme_void()
a  <- plot_grid(hist_x, p, ncol=1, rel_heights=c(1,6), align="v")
b  <- plot_grid(empty, hist_y, ncol=1, rel_heights=c(1,6))
p_final <- plot_grid(a, b, rel_widths=c(6,1))
p_final

#fig7d,7e
library(DESeq2)
library(Seurat)
library(IHW)
library(monocle)
library(tidyverse)
library(magrittr)
library(Seurat) 
library(pheatmap)
library(RColorBrewer)
library(aplot)

RA_matrix<-as(as.matrix(data.combined@assays$RNA@counts), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<- data.combined@meta.data
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
monocle_cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- setOrderingFilter(monocle_cds, degs_sig)
disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.2)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
monocle_cds <- reduceDimension(monocle_cds, 
                               reduction_method = "DDRTree",
                               num_dim = 5, 
                                scaling=T,
                               auto_param_selection = T)
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "celltype", cell_size = 1.5)+
  scale_colour_manual(values = celltype_col)
ggsave("monocle2_cell.pdf", width=4, height=4)

plot_cell_trajectory(monocle_cds, color_by = "group", cell_size = 1.5)+
  scale_colour_manual(values = sample_color)
ggsave("monocle2_group.pdf", width=4, height=4)

plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1.5)+
  scale_color_gradientn(colors  = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
ggsave("monocle2_time.pdf", width=4, height=4)

plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1.5)+
  scale_colour_manual(values = col)
ggsave("monocle2_State.pdf", width=4, height=4)
branch_point = 1

BEAM_res=BEAM(monocle_cds[expressed_genes,],branch_point = branch_point,cores = 13)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")
tmp1=plot_genes_branched_heatmap(monocle_cds[row.names(BEAM_res)[1:60]],
                                 branch_point = branch_point,
                                 num_clusters = 2, 
                                 cores = 13,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, 
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("grey80", "#F05662", "#FF935F"), 
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap = T )

pdf("branched_heatmap.pdf",width = 4.8,height = 5.2)
tmp1$ph_res
dev.off()

#fig7f
degs <- FindMarkers(object = data.combined, 
                       only.pos = TRUE,
                       ident.1="B",
                       ident.2="A",
                       logfc.threshold = 0.1) 

diff=degs
diff$P.Value=diff$p_val
diff$logFC=diff$avg_log2FC
write.csv(diff,"diff.csv")
diff_sub=diff[diff$P.Value<0.05,]
diff_sub=diff_sub[diff_sub$logFC>0.05,]
degs_all=unique(rownames(diff_sub))
sigScores <- as.matrix(data.combined@assays$RNA@counts[degs_all,])
sigScores_degs=sigScores[degs_all,]
sigScores_degs[which(sigScores_degs>0)]<-1
degs_express=colSums(sigScores_degs)
meta=data.combined@meta.data
meta=meta[names(degs_express),]
meta$degs_express=degs_express
library("viridis")
UMap = data.frame(data.combined@reductions$umap@cell.embeddings)
UMap$TranscriptCount = meta$degs_express
UMap$Cluster = meta$fincell
max(UMap$TranscriptCount)
UMap$TranscriptCount[which(UMap$TranscriptCount>435)]<-435
ggplot(UMap[order(UMap$TranscriptCount),],aes(x=UMAP_1, y=UMAP_2, colour=TranscriptCount)) + 
  geom_point(alpha=0.8, size=0.9) +
  xlab("UMAP 1") + 
  ylab("UMAP 2") + 
  theme_bw()+
  theme(panel.border = element_rect(fill=NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        aspect.ratio=1,  # width : height =1 
        legend.position = 'right',
        legend.key = element_rect(size=1),# legend size
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
  ) + 
  labs(colour="Number of genes")+
  scale_color_viridis(option = "C", direction = 1)

#fig7g
library(ggplot2)
library(jjPlot)
meta.data=data.combined@meta.data
data=data.frame(meta.data$celltype,meta.data$origin,meta.data$group)
colnames(data)=c("cell","orin","group")
df.long <- reshape2::melt(table(data))
df.long_1=df.long[df.long$orin=="Before",]
df.long_2=df.long[df.long$orin=="Post",]
colnames(df.long_2)=c("cell2","orin2","time2","value2")
df=cbind(df.long_1,df.long_2)
df$id=rownames(df)
df$group=rep(c(1:2),c(7,7))
ggplot(df,aes(x = time,group = group,y=1)) +
  geom_jjPointPie(aes(pievar = value2,
                      fill = cell2),
                  width = 1.4,
                  color = "white",
                  add.text = T,
                  text.dist = 0.3,
                  text.size = 0) +
  geom_tile(fill = NA,color = 'black') +
  coord_fixed() +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

#fig7h
library(stringr)
data.combined$sam_group=paste0(data.combined$sam,":",data.combined$group)
cellratio <- prop.table(table(data.combined$fincell, data.combined$sam_group), margin = 2)%>%
  as.data.frame()%>%
  `colnames<-`(c('celltype','orig.ident','Freq'))
group=strsplit(as.vector(cellratio$orig.ident),":")
group=lapply(group,FUN = function(x){x[2]})
group=unlist(group)
cellratio$group <-group
id=strsplit(as.vector(cellratio$orig.ident),":")
id=lapply(id,FUN = function(x){x[1]})
id=unlist(id)
cellratio$id <-id
df=cellratio
ggplot(df, aes(x = group, y = Freq, color = group)) +
  geom_line(aes(group = id), color = "grey80", size = 0.5) +
  geom_point(size = 1) +
  stat_compare_means(paired = TRUE, size = 2) +
  facet_wrap(~ celltype, scales = 'free_y', nrow = 2) +
  scale_fill_manual(values=sample_color)+
  scale_color_manual(values=sample_color)+
  theme_classic2(base_size = 9) +
  facet_theme +
  xlab('Freq')
ggsave('pair_all4.pdf', width = 10, height = 4)

#fig7i
library(ggtern)
library(scales)
library(ggplot2)
library(limma)
sigScores <- as.matrix(t(getSignatureScores(vis)))
de_gsva <- function(exprSet,meta,compare = NULL,olny_two=T){
   allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  fit <- lmFit(exprSet,design)
  if (olny_two==T) { 
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))==21){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}
Diff =de_gsva(exprSet = sigScores ,meta = meta,compare = compare,olny_two=F)
Padj_threshold=0.01
logfc=0.06
path=c()
for (i in c(1:length(Diff))) {
  data=Diff[[i]]
  tmp=data[data$adj.P.Val<Padj_threshold,]
  tmp=tmp[abs(tmp$logFC)>logfc,]
  tmp=rownames(tmp)
  path=c(path,tmp)
}
path=unique(path)
library("Seurat")
library("pheatmap")
library("viridis")
library("dplyr")
avgData <- sigScores %>% 
  apply(1, function(x){
    tapply(x, meta, mean) # ExpMean
  }) %>% t
avgData2=avgData[path,]
phData=avgData2
head(phData)
data=phData[,c(1:3)]
library(ClassDiscovery)
hcs <- hclust(distanceMatrix(as.matrix(t(data)), 'sqrt pearson'), "ward.D")
scgroup <- cutree(hcs,k = 3) 
scgroup=as.factor(scgroup)
data=as.data.frame(data)
data$group=scgroup
ggtern(data=data, aes(x=A,
                      y=D,
                      z=B_C,color=group))+geom_point(aes(size=0.5), alpha=0.5)+ scale_color_manual(values =sample_color2)+
  theme(tern.panel.background = element_rect(fill = "white"),
        tern.panel.grid.minor = element_line(color = "gray90"), 
        tern.axis.arrow.show = TRUE, 
        tern.axis.arrow.T = element_line(color ='#FB7D80', size = 0.1), 
        tern.axis.arrow.L = element_line(color = 'grey30', size = 0.1),
        tern.axis.arrow.R = element_line(color = '#6FC7CF', size = 0.1),
        tern.axis.arrow.text.L = element_text(color = 'black'),  
        tern.axis.arrow.text.T = element_text(color = 'black'),
        tern.axis.arrow.text.R = element_text(color = 'black'),
        tern.axis.arrow.sep = 0.1, 
        tern.panel.grid.major.T = element_line(color = 'gray92', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.major.L = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.major.R = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.T = element_line(color = 'gray94', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.minor.L = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.R = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.axis.title.L = element_text(color = 'grey30', size = 11),
        tern.axis.title.T = element_text(color = '#FB7D80', size = 11),
        tern.axis.title.R = element_text(color = '#6FC7CF', size = 11),
        tern.axis.text.L = element_text(size = 17,face = 'bold'),
        tern.axis.text.R = element_text(size = 17,face = 'bold'),
        tern.axis.text.T = element_text(size = 17,face = 'bold'),
        tern.axis.vshift = 0.04,
        tern.axis.line.T = element_line(linewidth = 0.8),
        tern.axis.line.R = element_line(linewidth = 0.8),
        tern.axis.line.L = element_line(linewidth = 0.8))

