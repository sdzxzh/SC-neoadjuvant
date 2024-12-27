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

#fig5a
#refer to fig4b

#fig5b
#refer to fig4c

#fig5c
#refer to fig4d
library(LSD)
pdf(paste0("heatscatter","_",name,".pdf"),4,4.5)
heatscatter(data.combined@reductions$umap@cell.embeddings[,1],
            data.combined@reductions$umap@cell.embeddings[,2])
dev.off()

#fig5e
library(viridis)
library(ggdensity)
library(ggblanket)
library(ggsci)
sample_color <- c("#1CC5FE","#FB7D80","#FBA27D","#39398A","#AF1E1F","#FBA27D")
gene_all<-c("CXCL13")

KS_plot_density2 <- function(obj,
                            marker,
                            dim=c("TSNE","UMAP"),
                            size,
                            ncol=NULL
){
  require(ggplot2)
  require(ggrastr)
  require(Seurat)
  
  cold <- colorRampPalette(c('#f7fcf0','#00bbb1','#0f7ab0'))
  warm <- colorRampPalette(c('#ffffb2','#ec840e','#ca443d'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  if(dim=="TSNE"){
    
    xtitle = "tSNE1"
    ytitle = "tSNE2"
    
  }
  
  if(dim=="UMAP"){
    
    xtitle = "UMAP1"
    ytitle = "UMAP2"
  }
  
  
  if(length(marker)==1){
    
    plot <- FeaturePlot(obj, features = marker)
    data <- plot$data
    
    
    if(dim=="TSNE"){
      
      colnames(data)<- c("x","y","ident","gene")
      
    }
    
    if(dim=="UMAP"){
      
      colnames(data)<- c("x","y","ident","gene")
    }
    p <- ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,], 
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size,alpha=0.8) +
      geom_point_rast(data=data[data$gene>0,], 
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size,alpha=0.8) +
      scale_fill_gradientn(colours = mypalette)+
      scale_colour_gradientn(colours = mypalette)+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) 
    return(p)
  }else{
    gene_list <- list()
    for (i in 1:length(marker)) {
      plot <- FeaturePlot(obj, features = marker[i])
      data <- plot$data
      if(dim=="TSNE"){  
        colnames(data)<- c("x","y","ident","gene")
      }
      if(dim=="UMAP"){   
        colnames(data)<- c("x","y","ident","gene")
      }
      gene_list[[i]] <- data
      names(gene_list) <- marker[i]
    }
    plot_list <- list()
    for (i in 1:length(marker)) {
      p <- ggplot(gene_list[[i]], aes(x, y)) +
        geom_point_rast(shape = 21, stroke=0.25,
                        aes(colour=gene, 
                            fill=gene), size = size,alpha=0.8) +
        geom_density_2d(data=gene_list[[i]][gene_list[[i]]$gene>0,], 
                        aes(x=x, y=y), 
                        bins = 5, colour="black") +
        scale_fill_gradientn(colours = mypalette)+
        scale_colour_gradientn(colours = mypalette)+
        theme_bw()+ggtitle(marker[i])+
        labs(x=xtitle, y=ytitle)+
        theme(
          plot.title = element_text(size=12, face="bold.italic", hjust = 0),
          axis.text=element_text(size=8, colour = "black"),
          axis.title=element_text(size=12),
          legend.text = element_text(size =10),
          legend.title=element_blank(),
          aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      plot_list[[i]] <- p
    }
    Seurat::CombinePlots(plot_list, ncol = ncol)
  }
}
Idents(data.combined)<-"group"
obj=data.combined
gene_threshold=1
gene_threshold_density=1
gene_all=intersect(gene_all,rownames(data.combined@assays$RNA@counts))
for (i in gene_all) {
  
  if (sum(obj@assays$RNA@counts[i,])!=0) {

    dim="UMAP"
    marker=i
    cold <- colorRampPalette(c('#f7fcf0','#00bbb1','#0f7ab0'))
    warm <- colorRampPalette(c('#ffffb2','#ec840e','#ca443d'))
    mypalette <- c(rev(cold(11)), warm(10))
    
    if(dim=="UMAP"){
      
      xtitle = "UMAP1"
      ytitle = "UMAP2"
    }
    p<-KS_plot_density(obj=obj, 
                       marker= i,
                       dim = "UMAP", size =0.9)
    
    ggsave(plot=p,filename=paste0(i,"_type1.pdf"),width = 3.5, height = 3.5)
    p<-KS_plot_density2(obj=obj, 
                       marker= i,
                       dim = "UMAP", size =0.9)
    
    ggsave(plot=p,filename=paste0(i,"_type1_2.pdf"),width = 3.5, height = 3.5)
    plot <- FeaturePlot(obj, features = i)
    data <- plot$data
    colnames(data)<- c("x","y","ident","gene")
    p2<-ggplot(data[data$gene>gene_threshold_density,],aes(x = x, y = y, fill = ident)) +
      geom_hdr( )+
      xlim(min(p$data$x),max(p$data$x))+
      ylim(min(p$data$y),max(p$data$y))+
      #scale_fill_aaas() +
      #scale_fill_discrete(type  = c('#ffffb2','#ffffb2','#ec840e','#ca443d'))+
      scale_fill_manual(values = sample_color)+
      theme_bw()+
      facet_wrap(vars(ident))+
      labs(x="UMAP_1", y="UMAP_2")+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    ggsave(plot=p2,filename=paste0(i,"_type2.pdf"),width = 8, height = 4)
    
    p3<-ggplot(data[data$gene>gene_threshold_density,],aes(x = x, y = y)) +
      geom_hdr( )+
      xlim(min(p$data$x),max(p$data$x))+
      ylim(min(p$data$y),max(p$data$y))+
      theme_bw()+
      labs(x="UMAP_1", y="UMAP_2")+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    ggsave(plot=p3,filename=paste0(i,"_type3.pdf"),width = 4, height = 4)
    p4<-FeaturePlot(obj, raster=F,
                    
                    #cols=c('grey',"#FFFF99","#CC3333"),
                    # cols=c("darkblue", "lightblue", "green", "yellow", "red"),
                    features =c(i)
                    
    )
    data=p4$data
    colnames(data)<- c("x","y","ident","gene")
    p4<-ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      geom_point_rast(data=data[data$gene>gene_threshold,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      
      scale_fill_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                       '#bef0b0', '#fdf4af', '#f9b64b',
                                       '#ec840e', '#ca443d', '#a51a49'))+
      scale_colour_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                         '#bef0b0', '#fdf4af', '#f9b64b',
                                         '#ec840e', '#ca443d', '#a51a49'))+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    ggsave(plot=p4,filename=paste0(i,"_type4_1.pdf"),width = 3.5, height = 3.5)
    
    p4<-p4+geom_density_2d(data=data[data$gene>gene_threshold,], 
                           aes(x=x, y=y), 
                           bins = 5, colour="black") 
    ggsave(plot=p4,filename=paste0(i,"_type4_2.pdf"),width = 3.5, height = 3.5)
    p5<-ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      geom_point_rast(data=data[data$gene>gene_threshold,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    p5<-p5+scale_color_viridis( option = "D")+scale_fill_viridis( option = "D")
    
    ggsave(plot=p5,filename=paste0(i,"_type5_1.pdf"),width = 3.5, height = 3.5)
    
    p5<-p5+geom_density_2d(data=data[data$gene>gene_threshold,], 
                           aes(x=x, y=y), 
                           bins = 5, colour="black") 
      ggsave(plot=p5,filename=paste0(i,"_type5_2.pdf"),width = 3.5, height = 3.5)
    p6<-ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      geom_point_rast(data=data[data$gene>gene_threshold,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    p6<-p6+scale_color_viridis( option = "A")+scale_fill_viridis( option = "A")
    ggsave(plot=p6,filename=paste0(i,"_type6_1.pdf"),width = 3.5, height = 3.5)
    
    p6<-p6+geom_density_2d(data=data[data$gene>gene_threshold,], 
                           aes(x=x, y=y), 
                           bins = 5, colour="black") 
      ggsave(plot=p5,filename=paste0(i,"_type6_2.pdf"),width = 3.5, height = 3.5)
    p7<-ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      geom_point_rast(data=data[data$gene>gene_threshold,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    p7<-p7+scale_color_viridis( option = "F")+scale_fill_viridis( option = "F")
    ggsave(plot=p7,filename=paste0(i,"_type7_1.pdf"),width = 3.5, height = 3.5)
    p7<-p7+geom_density_2d(data=data[data$gene>gene_threshold,], 
                           aes(x=x, y=y), 
                           bins = 5, colour="black") 
    ggsave(plot=p7,filename=paste0(i,"_type7_2.pdf"),width = 3.5, height = 3.5)
    p8<-ggplot(data, aes(x, y)) +
      geom_point_rast(data=data[data$gene==0,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      geom_point_rast(data=data[data$gene>gene_threshold,],
                      shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = 0.9,alpha=0.8) +
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    p8<-p8+scale_color_gradient(low = "lightgrey", high = "#D24641")+scale_fill_gradient( low = "lightgrey", high = "#D24641")
    ggsave(plot=p8,filename=paste0(i,"_type8_1.pdf"),width = 3.5, height = 3.5)
    p8<-p8+geom_density_2d(data=data[data$gene>gene_threshold,], 
                           aes(x=x, y=y), 
                           bins = 5, colour="black") 
    ggsave(plot=p8,filename=paste0(i,"_type8_2.pdf"),width = 3.5, height = 3.5)
    
    
  }
  
}

#fig6k,l
library(survminer)
library(survival)
risk=dat$gene
res.cut=surv_cutpoint(dat,time="OS.time",
                      event ="OS",variables="sig")
res.cut=res.cut$cutpoint$cutpoint

risk<-as.vector(ifelse(risk >res.cut,"high","low"))
dat$group<-risk
sur=dat
fitd <- survdiff(Surv(OS.time,OS) ~ group,
                 data      = dat,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = dat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

ps <- pairwise_survdiff(Surv(OS.time, OS)~ group,
                        data = dat) 
names(fit$strata) <- gsub("group=", "", names(fit$strata))

p <- ggsurvplot(fit = fit,
                conf.int = F, risk.table= T, # 生存风险表risk.table.col    = "strata",
                palette = mycol, data= dat,size= 1,
                break.time.by= 5,legend.title= "",
                xlab= "Time (years)",
                ylab= "Overall survival",
                risk.table.y.text = FALSE,
                tables.height     = 0.3) 
p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", 
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55, 
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
p

#figS6a
library(ggplot2)
p <- DotPlot(data.combined_fib, features = rev(gene),#features = top5$gene,
             cols = my36colors, 
             assay = "RNA",
             group.by = "fincell", split.by = "group")+coord_flip()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)

exp$id <- factor(exp$id,levels = levels(p$data$id))
p1 <- ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="white",
             stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=11,color="black"),
        panel.margin=unit(c(2,0,0,1), "cm"),
        plot.margin = margin(t = 5,r = 1,b = 1,l = 1,unit = 'cm'),
        axis.text.x = element_text(angle=45,
                                   vjust=1, 
                                   size=11,
                                   hjust=1,
                                   
                                   color = "black")
  )+
  coord_cartesian(clip = 'off') +
  scale_color_gradientn(colors = colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)+ 
  geom_vline(xintercept=c(3.5,6.5,9.5), linetype="dotted",size=1.2)
p2 <- annoSegment(object = p1,
                  annoPos = 'top',
                  annoManual = F,
                  xPosition = c(1:13),
                  yPosition = c(59.5),
                  segWidth = 0.7,
                  pCol=rep(c(rep("#57C3F3",1),rep("#FB7D80",1)),3)
)

p3 <- annoSegment(object = p2,
                  annoPos = 'top',
                  annoManual = T,
                  xPosition = list(c(1,3,5,7),
                                   c(2,4,6,8)),
                  yPosition = c(48.5),
                  segWidth = 0.7,
                  pCol = col,
                  addBranch = T,
                  branDirection = -1,
                  lwd = 3)

#figS6b
library(ggplot2)
library(dplyr)
library(reshape2)
bar_plot_input$N=-bar_plot_input$N
x=melt(bar_plot_input,id.vars = "cell")
x$ylab=x$cell
x$ylab=factor(x$ylab,levels = bar_plot_input[order(bar_plot_input$T),]$cell)
x$variable=factor(x$variable,levels = c("T","N"))
x$variable=factor(x$variable,levels = c("T","N"))
ggplot(x,aes(x=ylab,y=value,fill=variable))+
  geom_col(position=position_dodge(0),width = 1)+
  coord_flip()+
  ylim(-1.65,1.65)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  scale_fill_manual(guide = guide_legend(title = NULL),values = col)+
  theme(legend.position ="top", 
       # axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(), 
        legend.key.size=unit(0.8,"line"),
        axis.text.y = element_text(color =celltype_col[levels(x$ylab)],size = 12))
  
#figS6c
#refer to fig5f









