library(tidyr)
library(dplyr)
library(tibble)
library(ggrepel)
library(tidyverse)
library(limma)
library(edgeR)
library(pheatmap)
setwd("/home/leili/GSEBJP")
#####------------------------C.heatmap----------------------------
read1 <- read_delim("GSE217853_Ctrl-fpkm.txt", delim = "\t")
read2 <- read_delim("GSE217853_DM-fpkm.txt", delim = "\t")
read <- left_join(read1, read2, by = "id") %>%
  column_to_rownames("id")
logread <- edgeR::cpm(read, log = TRUE)
heatmap <- logread[geneList0,]
treat <- data.frame(sample = colnames(logread), gt = rep(c("Ctrl", "DM"), each = 18),
                    time = rep(c("ZT0", 'ZT4', 'ZT8', 'ZT12', 'ZT16', 'ZT20'), each =3))
gt <- factor(treat$gt)

anno_df <- treat %>% 
  column_to_rownames('sample') %>% 
  rename(group = gt)
pheatmap(heatmap, cluster_cols = FALSE,cluster_rows = FALSE,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         scale = 'row', cellheight = 10, 
         annotation_col = anno_df[,1, drop = FALSE])

#####---------------------------E.volcano plot--------------------------
result <- read.csv("/home/leili/GSEBJP/result.csv",row.names = 1)
result <- rename(result, log2FoldChange = logFC, padj = adj.P.Val)

result$change = ifelse(result$padj < 0.05 & abs(result$log2FoldChange) >= 0, 
                       ifelse(result$log2FoldChange> 0 ,'Up','Down'),
                       'Stable')

#需要突出显示的基因列表(用户自定义)
geneList0 <- c('Adrb1',"Adrb2",'Adrb3',
               'Pde1a','Pde4d','Pde4b','Pde3a','Pde5a','Pde8a',
               'Pde4a','Pde1b','Pde1c','Pde6d',
               'Ccn2', 'Acta2', "Maoa",
               'Tgfbr3l','Tgfb1i1','Tgfbr1')
geneList <- result[geneList0,]

library('ggplot2')
p <- ggplot(# 数据、映射、颜色
  result, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #突出表示差异基因
  geom_point(data=geneList,aes(x = log2FoldChange, y = -log10(padj)),colour="yellow",size=1.5)+
  #辅助线
  geom_vline(xintercept=c(-0,0),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))
p

#标记出5个基因的label
geneList1 <- result[rownames(result) %in% geneList0,]
geneList1 <- subset(geneList1, select = -change)
geneList1$label <- rownames(geneList1)

library(ggrepel)
p + geom_label_repel(data = geneList1, 
                     aes(x = log2FoldChange, y = -log10(padj), label = label),
                     size = 3,color="black",
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",   #连线的颜色
                     segment.size = 0.4,max.overlaps = 30)

#####----------------------------F.Bpxplot ADRB1---------------------------
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
read1 <- read_delim("GSE217853_Ctrl-fpkm.txt", delim = "\t")
read2 <- read_delim("GSE217853_DM-fpkm.txt", delim = "\t")
read <- left_join(read1, read2, by = "id") %>%
  column_to_rownames("id")
#logread <- edgeR::cpm(read, log = TRUE)

Exp <- log2(read+1)
anno <- data.frame(sample = colnames(Exp),group = rep(c("Ctrl","DM"),each = 18),
                   time = rep(c("ZT0", 'ZT4', 'ZT8', 'ZT12', 'ZT16', 'ZT20'), each =3))
Exp <- t(Exp)
#anno <- anno %>% 
#  column_to_rownames('sample')

gene <- c("Adrb1",
          'Pde4d','Pde4b')#这里我们只选择这几个基因做数据
gene <- as.vector(gene)
Exp_plot <- Exp[,gene]#提取需要作图得基因表达信息
#加载样本信息
Exp_plot <- Exp_plot[anno$sample,]
Exp_plot <- as.data.frame(Exp_plot)
Exp_plot$sam = anno$group
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("Ctrl","DM"))
######
library("dplyr")
group_by(Exp_plot, sam) %>%
  summarise(
    count = n(),
    mean = mean(Adrb1, na.rm = TRUE),
    sd = sd(Adrb1, na.rm = TRUE)
  )
ggboxplot(Exp_plot, x = "sam", y = "Adrb1", 
          color = "sam", palette = c("#00AFBB", "#E7B800"),
          order = c("Ctrl","DM"),
          ylab = "Adrb1", xlab = "Groups")
ggboxplot(Exp_plot, x = "sam", y = "Pde4d", 
          color = "sam", palette = c("#00AFBB", "#E7B800"),
          order = c("Ctrl","DM"),
          ylab = "Pde4d", xlab = "Groups")
#####
col <-c("#337AB7","#D9534F","#5CB85C")#"#5CB85C","#F0AD4E",
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp <- Exp_plot[,c(gene[i],"sam")]
  colnames(bar_tmp) <- c("Expression","sam")
  my_comparisons <- list(c("Ctrl","DM"))
  pb1 <- ggboxplot(bar_tmp,
                   x="sam",
                   y="Expression",
                   color="sam",
                   fill=NULL,
                   add = "jitter",
                   bxp.errorbar.width = 0.6,#0.6
                   width = 0.4,
                   size=0.01,
                   font.label = list(size=30), 
                   palette = col)+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")#
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = T,
                              comparisons =c(my_comparisons),
                              label="p.format")
  plist2[[i]]<-pb1 
}
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],ncol=4)#ncol=4表示图片排为几列
