rm(list = ls())
options(stringsAsFactors = FALSE)

# 加载R包
library(Seurat)
library(harmony)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggsci)
library(dplyr)
library(tidyverse)
library(magrittr)
library(scRepertoire)
library(coin)
library(patchwork)
library(stringr)
library(cowplot)
library(DESeq2)
library(scales)
library(ggpubr)

major_colors <- c("CD4T cells" = "#4f6980",
                  "CD8T cells" = "#8397a6",
                  "NK cells" = "#bfbb60",
                  "B&Plasma Cell" = "#a2ceaa",
                  "Mast cell" = "#b66353",
                  "Myeloid" = "#e0c8bd",
                  "Epithelia Cell" = "#fbb04e",
                  "Fibroblast" = "#f47942")
# fig1b
plot <- DimPlot(seu, reduction = "umap", group.by = "cellType", label = F, raster = FALSE, cols = major_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL) 
plot
ggsave(
  filename = "seu_dimplot.eps",   
  plot = plot,          
  device = "eps",             
  width = 8,                  
  height = 6,                 
  units = "in",             
  dpi = 300            
)

sample_colors <- c("CRN" = "#3d7dba",
                   "CRT" = "#c6a5cb",
                   "PLN" = "#3ba570",
                   "WBC" = "#f7c497")
plot <- DimPlot(seu, reduction = "umap", group.by = "sampleType", label = F, raster = FALSE, cols = sample_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)
plot
ggsave(
  filename = "seu_sample_dimplot.eps",   
  plot = plot,       
  device = "eps",            
  width = 8,              
  height = 6,                
  units = "in",            
  dpi = 300     
)

# fig1c
all_markers <- unique(c("PTPRC","CD3E","CD8A","CD4","TYROBP","TRGV9","TRDV2",#CD8T_markers,CD4T_markers,
                        "CD19","MZB1",#Bpls_marker,
                        "CD14","CD68","FCGR3A",#macrophage_markers,
                        "CLEC9A","XCR1","BATF3","CLEC10A","CD1C","LILRA4","IL3RA","CLEC4C",
                        "TPSAB1",
                        "EPCAM","VWF","COL1A1"))
DotPlot(seu, features = all_markers, group.by = "cellType") + coord_flip()
# 使用reorder函数来调整横坐标的顺序
reordered_annotations <- c("CD4T cells","CD8T cells","NK cells","B&Plasma Cell","Myeloid","Mast cell","Fibroblast","Epithelia Cell")
# 重新排列横坐标的顺序
seu$cellType <- factor(seu$cellType, levels = reordered_annotations)
# 创建DotPlot
DotPlot(seu, features = all_markers, group.by = "cellType") +
  # coord_flip() +
  labs(
    title = "",  # 设置图的标题
    x = "",     # 设置x轴标签
    y = ""      # 设置y轴标签
  ) +
  scale_color_gradient2(low = "#4575b4", high = "#d73027") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# fig1d
imm$cellType <- factor(imm$cellType, levels = reordered_annotations)
tmp_info <- select(imm@meta.data, "sampleType", "cellType")
tmp_info$sampleType <- factor(tmp_info$sampleType, levels = c("PLN", "WBC", "CRN", "CRT"))
tmp_info %>%
  dplyr::count(sampleType, cellType) %>% 
  group_by(sampleType) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=sampleType, y=proportion, fill=cellType)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Sample Type", fill="Cell Class") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = major_colors)

# fig1e
cd8t_colors <- c("CD8T_CCR7" = "#87c097",
                 "CD8T_GZMK" = "#166a3b",
                 "CD8T_FGFBP2" = "#00b8c4",
                 "CD8T_CXCL13" = "#c571ac",
                 "CD8T_HSPA1A" = "#c2aed1",
                 "CD8T_CD160" = "#9775b7",
                 "CD8T_MAIT" = "#feac75")

DimPlot(cd8t_1, reduction = "umap", group.by = "cell_annotation", label = F, raster = FALSE, cols = cd8t_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)

# fig1f
Idents(cd8t_1) <- cd8t_1$cell_annotation
CD8T_markers <- c("SELL","LEF1","TCF7","CCR7",
                  "LAG3","TIGIT","PDCD1","CTLA4","ENTPD1","LAYN","HAVCR2",
                  "PRF1","GZMB","GZMK","NKG7","GZMA","KLRD1","FCRL6","GNLY","IFNG","TNF","FGFBP2",
                  "ICOS","CD28","CD27","TNFRSF9","TNFRSF18",
                  "CXCL13","BHLHE40","ITGAE","CD69","CD160",
                  "HSPA1A","HSPA1B","HSP90AA1","HSP90AB1","ISG15")
cd8t_matrix <- AverageExpression(cd8t_1, slot = "scale.data")[[1]]  
order_index <- match(CD8T_markers, rownames(cd8t_matrix))
col_index <- match(c("CD8T_CCR7","CD8T_GZMK","CD8T_FGFBP2","CD8T_CD160","CD8T_MAIT","CD8T_HSPA1A","CD8T_CXCL13"),colnames(cd8t_matrix))
cd8t_matrix <- cd8t_matrix[order_index, ]
cd8t_matrix <- cd8t_matrix[, col_index]
cd8t_matrix[cd8t_matrix > 2] = 2;cd8t_matrix[cd8t_matrix < -2] = -2

pheatmap::pheatmap(cd8t_matrix, scale = "row", cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7", "#F3DBD6", "#DA8F87", "#D54846"))(100),
                   gaps_row = c(4,11,22,27),
                   gaps_col = c(2,3,4,5),
                   # annotation_row = annotation.row.data,
                   # annotation_colors = annotation_colors
) #4,9

cd8t_info <- cd8t_1@meta.data
clusters_SID_ob <- as.matrix(table(cd8t_info$sampleType, cd8t_info$cell_annotation))
cellsum <- table(cd8t_info$cell_annotation)
SIDsum <- table(cd8t_info$sampleType)
clusters_SID_exp <- matrix(rep(cellsum, length(SIDsum)), nrow = length(SIDsum), byrow = TRUE) * matrix(rep(SIDsum, length(cellsum)), nrow = length(SIDsum), byrow = FALSE) / nrow(cd8t_info)
Ro_e <- clusters_SID_ob/clusters_SID_exp
# range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
# Ro_e_01 <- t(apply(Ro_e, 1, range01))
bk <- c(seq(0, 1, by = 0.01), seq(1.01, 2, by = 0.01))
col_index <- match(c("CD8T_CCR7","CD8T_GZMK","CD8T_FGFBP2","CD8T_CD160","CD8T_MAIT","CD8T_HSPA1A","CD8T_CXCL13"),colnames(Ro_e))
Ro_e <- Ro_e[, col_index]

Ro_e %>%
  pheatmap::pheatmap(cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     scale = "none",
                     color = c(colorRampPalette(colors = c("#FFFFB2","#FD8D3C"))(length(bk)/2),
                               colorRampPalette(colors = c("#FD8D3C","#B10026"))(length(bk)/2)),
                     breaks=bk,
                     border_color = NA, 
                     display_numbers = T,
                     annotation_legend = FALSE,
                     number_color = "grey20", 
                     fontsize_number = 10,
                     gaps_col = c(2,3,4,5), show_colnames = F)

# fig1g
tmp_info <- select(cd8t_1@meta.data, "sampleType", "cell_annotation")
# 设置 sampleType 的顺序
tmp_info$sampleType <- factor(tmp_info$sampleType, levels = c("PLN", "WBC", "CRN", "CRT"))
tmp_info$cell_annotation <- factor(tmp_info$cell_annotation, levels = rev(c("CD8T_CCR7","CD8T_GZMK","CD8T_FGFBP2","CD8T_CD160","CD8T_MAIT","CD8T_HSPA1A","CD8T_CXCL13")))
tmp_info %>%
  dplyr::count(sampleType, cell_annotation) %>% 
  group_by(sampleType) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=sampleType, y=proportion, fill=cell_annotation)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Sample Type", fill="cell annotation") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = cd8t_colors)

tmp_info <- select(cd8t_1@meta.data, "sampleType", "exhaustion.score")
tmp_info$sampleType <- factor(tmp_info$sampleType, levels = c("WBC", "PLN", "CRN", "CRT"))
compaired <- list(c("CRN","CRT"),
                  c("PLN","CRT"),
                  c("WBC","CRT"))
ggplot(data = tmp_info, aes(x = sampleType, y = exhaustion.score, fill = sampleType)) +
  geom_violin(scale = "width", width = 0.8, alpha = 0.8) +
  # geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white") +
  # geom_point() +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Exhaustion score")+
  scale_x_discrete(name = "") +
  # labs(fill = "sample_type") +
  # ggtitle("Boxplot of CD8+T") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 10),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10, 
                                 angle = 45,
                                 hjust =1)) +
  scale_fill_manual(values = sample_colors)+
  geom_signif(comparisons = compaired,
              step_increase = 0.3,
              map_signif_level = T,
              test = wilcox.test)

# fig2a-c
ggplot(data = cd8t_tmp.df, aes(x = naive.score)) + 
  geom_density() +  # 不填充颜色
  # scale_color_manual(values = c(c("MSI" = "#ED0000FF",
  #                                 "MSS" = "#00468BFF"))) +  # 手动设置颜色
  labs(title = "Density Plot", x = "naive.score", y = "Density") +
  theme_classic()

ggplot(data = cd8t_tmp.df[cd8t_tmp.df$cell_class %in% c("Tnaive","Ttsm","Tpex"),], aes(x = naive.score, y = prog.score)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.8, bins = 7, contour = TRUE) +  # 使用after_stat修正，调整透明度
  scale_fill_gradientn(colors = c("#fcf1f0", "#fccccb", "#e3716e"),  # 设置颜色渐变
                       values = scales::rescale(c(0, 0.5, 1)),  # 调整颜色在密度级别上的位置
                       limits = c(0, NA)) +  # 设置下限，使低密度区域有颜色
  geom_hline(yintercept = 0.14, linetype = "dashed", color = "black", size = 0.5) +
  geom_segment(aes(y = 0.14, yend = max(cd8t_tmp.df$prog.score), x = 0.3, xend = 0.3), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  scale_y_continuous(name = "Memory score") +
  scale_x_continuous(name = "Naive score")

ggplot(data = cd8t_tmp.df[cd8t_tmp.df$cell_class %in% c("Tex_int","Tex_term","Tex_eff"),], 
       aes(x = exhaustion_score, y = effect_score)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, bins = 6, alpha = 0.8) +  # 设置 bins 为 3
  scale_fill_gradientn(colors = c("#fcf1f0", "#fccccb", "#e3716e"),  # 设置颜色渐变
                       values = scales::rescale(c(0, 0.5, 1))) +  # 调整颜色在密度级别上的位置
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "black", size = 0.5) +
  geom_segment(aes(x = 0.9, xend = max(cd8t_tmp.df$exhaustion_score), y = 0.9, yend = 0.9), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  scale_y_continuous(name = "Effector score") +
  scale_x_continuous(name = "Exhaustion score")

# fig2d
Idents(cd8t_tmp) <- cd8t_tmp$cell_class
CD8T_markers <- c("SELL","LEF1","TCF7","CCR7","IL7R",
                  "GZMK","PRF1","GZMA","GZMB","NKG7","KLRD1","GNLY","IFNG","TNF","FGFBP2",
                  "LAG3","TIGIT","PDCD1","CTLA4","ENTPD1","LAYN","HAVCR2","TOX",
                  "ITGAE","CD69","CXCL13","BHLHE40","TNFRSF9","TNFRSF18","ISG15")
cd8t_matrix <- AverageExpression(cd8t_tmp, slot = "scale.data")[[1]]  
order_index <- match(CD8T_markers, rownames(cd8t_matrix))
cd8t_matrix <- cd8t_matrix[order_index, ]
col_index <- match(c("Tnaive","Ttsm","Tpex","Tex_int","Tex_term","Tex_eff"),colnames(cd8t_matrix))
cd8t_matrix <- cd8t_matrix[, col_index]
cd8t_matrix[cd8t_matrix > 2] = 2;cd8t_matrix[cd8t_matrix < -2] = -2

pheatmap::pheatmap(cd8t_matrix, scale = "row", cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7", "#F3DBD6", "#DA8F87", "#D54846"))(100),
                   gaps_row = c(5,15,23,25),
                   gaps_col = c(1,2,3,4,5),
                   border_color = "white"
                   # annotation_row = annotation.row.data,
                   # annotation_colors = annotation_colors
) #4,9

# fig2e
head(cd8t_tmp.df)
compaired <- list(c("Tex_int", "Tex_term"),
                  c("Tex_eff", "Tex_term"),
                  c("Tex_int", "Tex_eff"))
cd8t_tmp.df$cell_class <- factor(cd8t_tmp.df$cell_class, levels = c("Tnaive", "Ttsm", "Tpex", "Tex_int", "Tex_term", "Tex_eff"))
p1 <- ggplot(data = cd8t_tmp.df, aes(x = cell_class, y = effect_score, fill = cell_class)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    axis.text.x = element_blank(),  # 隐藏横轴文字
    axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ex_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = TRUE,
    test = wilcox.test
  )

p2 <- ggplot(data = cd8t_tmp.df, aes(x = cell_class, y = exhaustion_score, fill = cell_class)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    axis.text.x = element_blank(),  # 隐藏横轴文字
    axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ex_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = TRUE,
    test = wilcox.test
  )
combined_plot <- wrap_plots(p1,p2,ncol = 1,nrow = 2)

# fig2f
DimPlot(cd8t_tmp, reduction = "umap", group.by = "cell_class", label = F, raster = FALSE, cols = ex_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)

# fig2g
tmp_info <- select(cd8t_tmp@meta.data, "sampleType", "cell_class")
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = rev(c("Tnaive","Ttsm","Tpex","Tex_int","Tex_eff","Tex_term")))
tmp_info$sampleType <- factor(tmp_info$sampleType, levels = c("PLN", "WBC", "CRN", "CRT"))
tmp_info %>%
  dplyr::count(sampleType, cell_class) %>% 
  group_by(sampleType) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=sampleType, y=proportion, fill=cell_class)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Sample Type", fill="Cell Class") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = ex_colors)

# fig2j
library(monocle)
cd8t_tmp.1 <- cd8t_tmp
#表达矩阵
data <- as(as.matrix(cd8t_tmp.1@assays$RNA@counts), 'sparseMatrix')
#表型信息
pd <- new('AnnotatedDataFrame', data = cd8t_tmp.1@meta.data)
#基因信息
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#构建CDS对象
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
HSMM <- monocle_cds

#使用monocle2选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 0.5 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "cell_class") +
  facet_wrap(~ cell_class) +
  scale_color_manual(values = ex_colors)

# 提取细胞名称和拟时序值
cell_names <- rownames(pData(HSMM))  # 获取细胞名称
pseudotime <- pData(HSMM)$Pseudotime  # 获取拟时序值

# 创建一个数据框，包含细胞名称和拟时序值
cell_pseudotime <- data.frame(row.names = cell_names, Pseudotime = pseudotime)
cd8t_tmp <- AddMetaData(cd8t_tmp, metadata = cell_pseudotime)
head(cd8t_tmp@meta.data)

# 绘制基于 pseudotime 的 UMAP 图
FeaturePlot(cd8t_tmp, features = "Pseudotime", reduction = "umap") +
  scale_color_gradient(low = "#000033", high = "#6699FF") +  # 从深蓝到浅蓝的颜色渐变
  guides(color = guide_colorbar(title = "Pseudotime")) +     # 添加颜色条标题
  labs(title = NULL) +                                       # 移除标题
  theme_classic()                                            # 使用简洁主题

# fig3a
# Tex-eff高表达的基因
table(cd8t_tmp$cell_class)
Idents(cd8t_tmp) <- cd8t_tmp$cell_class
tex_eff.deg <- FindMarkers(cd8t_tmp,
                           ident.1 = "Tex_eff", ident.2 = "Tnaive",
                           min.pct = 0.1,  
                           logfc.threshold = 0.1)
tex_eff.deg.filter <- tex_eff.deg %>% filter(p_val_adj < 0.05)
tex_eff.deg.filter <- tex_eff.deg.filter[order(tex_eff.deg.filter$avg_log2FC, decreasing = T),]

tex_eff.deg.all <- FindMarkers(cd8t_tmp,
                               ident.1 = "Tex_eff",
                               min.pct = 0.25,  
                               logfc.threshold = 0.1)
tex_eff.deg.all.filter <- tex_eff.deg.all %>% filter(p_val_adj < 0.05)
tex_eff.deg.all.filter <- tex_eff.deg.all.filter[order(tex_eff.deg.all.filter$avg_log2FC, decreasing = T),]

tex_eff.deg.tmp <- rownames(tex_eff.deg.filter)[1:100]
tex_eff.deg.all.tmp <- rownames(tex_eff.deg.all.filter)[1:100]
tex_eff.gene <- intersect(tex_eff.deg.tmp,tex_eff.deg.all.tmp)
save(tex_eff.gene, file = "tex_eff.gene.RData")

library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggfun)

tex_eff.deg.all <- arrange(tex_eff.deg.all, desc(avg_log2FC))
tex_eff.deg.all.select <- tex_eff.deg.all

deg_result <- data.frame("rank" = 1:nrow(tex_eff.deg.all.select),
                         "log2fc" = tex_eff.deg.all.select$avg_log2FC,
                         "pvalue" = tex_eff.deg.all.select$p_val_adj,
                         "gene" = rownames(tex_eff.deg.all.select))

deg_result$Symbol <- ifelse(deg_result$gene %in% genes, deg_result$gene, NA)

library(ggplot2)
library(ggrepel)
library(dplyr)

p <- deg_result %>%
  ggplot() + 
  geom_point(aes(x = rank, y = log2fc, color = pvalue, size = abs(log2fc)), alpha = 0.8) + 
  scale_color_gradientn(
    colors = c("#3d7dba","#a6c3db","#ffffff","#f3b8ac","#fb8072"),  # 自定义蓝-白-红渐变
    values = scales::rescale(c(1,0.75,0.5,0.25,0)),     # 调整渐变分布，使 p < 0.05 更偏向红色
    name = "P-value"
  ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
  geom_text_repel(data = deg_result,
                  aes(x = rank, y = log2fc, label = Symbol),
                  size = 2.5,
                  box.padding = 0.3,
                  # point.padding = 0.5,
                  segment.color = "grey",
                  segment.size = 0.1,
                  max.overlaps = 50,
                  # segment.angle = 30,    # 增加线段的角度
                  # segment.curvature = -0.1, # 增加线段的弯曲度
                  force = 1,              # 增强推开力
                  # force_pull = 1,          # 吸引力
                  direction = "both"
  ) + 
  scale_size(name = "log2(Fold Change)") + 
  labs(x = "Rank of Differentially Expressed Genes",
       y = "log2(Fold Change)") + 
  theme_bw() + 
  theme(
    # panel.border = element_rect(linewidth = 1, color = "black"),
    legend.background = element_rect(color = "white", linetype = 1),
    axis.text = element_text(color = "#000000", size = 8),
    axis.title = element_text(color = "#000000", size = 10)
  )
p

# fig3b
library(clusterProfiler)
library(org.Hs.eg.db)
ego_BP <- enrichGO(gene          = tex_eff.gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
ego_bp <- separate(data = ego_bp, 
                   col = GeneRatio,
                   into = c("GR1", "GR2"), 
                   sep = "/")
ego_bp <- mutate(ego_bp, 
                 GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
ego_bp <- ego_bp[order(ego_bp$GeneRatio, decreasing = TRUE),]
# ego_bp$Description <- factor(ego_bp$Description, levels = ego_bp$Description)

ego_bp_top10 <- ego_bp[1:10,]
ego_bp_top10 <- ego_bp_top10[order(ego_bp_top10$GeneRatio, decreasing = FALSE),]
ego_bp_top10$Description <- factor(ego_bp_top10$Description, levels = ego_bp_top10$Description)

ego_bp_top10$Description_wrapped <- str_wrap(ego_bp_top10$Description, width = 30)
# 绘图
ggplot(ego_bp_top10, aes(x = GeneRatio, y = Description)) +
  geom_segment(aes(x = 0, xend = GeneRatio, y = Description, 
                   yend = Description, color = -log10(p.adjust)), size = 1.5) + # 调整线段粗细
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_y_discrete(labels = ego_bp_top10$Description_wrapped) + # 使用换行后的标签
  scale_color_gradient(low = "#fdc5b8", high = "#d85152") +
  labs(
    x = "GeneRatio",
    y = NULL,
    # title = "Lollipop Plot for Pathway Enrichment",
    color = "-log10(P.adjust)",
    size = "Count"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10)
  )

# fig3c
tf.df <- rssPlot$df %>% dplyr::filter(Topic %in% c("ZNF121(+)","BATF(+)"))
ggplot(tf.df, aes(x = cellType, y = Topic, size = RSS, color = Z)) +
  geom_point() +
  scale_color_gradient2(
    low = "#4575b4", 
    mid = "white", 
    high = "#d73027", 
    midpoint = 1  # Set midpoint to the average Z value
  ) +
  scale_x_discrete(
    limits = c("Tex_eff", "Tex_term", "Tex_int", "Tpex", "Ttsm", "Tnaive")  # Set the order of x-axis
  ) +
  labs(size = "regulon specificity score", color = "Z score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  )

# fig3d
cd8texeff_info$Sample <- factor(cd8texeff_info$Sample, levels = c("Tumor","Invasion","Normal"))
compaired <- list(c("Tumor", "Invasion"),
                  c("Normal","Tumor"),
                  c("Invasion", "Normal"))
ggplot(data = cd8texeff_info, aes(x = Sample, y = Ratio_CD8, color = Sample)) +
  geom_boxplot(aes(fill = Sample), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Ratio of CD8_Texeff/CD8+T") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = c("Normal" = "#3d7dba",
                                "Tumor" = "#c6a5cb")) +
  scale_fill_manual(values = c("Normal" = "#3d7dba",
                               "Tumor" = "#c6a5cb")) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

ggplot(data = cd8texeff_info, aes(x = Sample, y = Density, color = Sample)) +
  geom_boxplot(aes(fill = Sample), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Density of CD8_Texeff cells/ROI") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = c("Normal" = "#3d7dba",
                                "Tumor" = "#c6a5cb")) +
  scale_fill_manual(values = c("Normal" = "#3d7dba",
                               "Tumor" = "#c6a5cb")) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig3g
flow_res.tmp <- dplyr::filter(flow_res, Cell_Group %in% c("TCF1-PD1hiGZMB+", "TCF1-PD1hiGZMB-"))
df_long <- flow_res.tmp %>%
  pivot_longer(
    cols = c("IL-2_Ratio", "TNF-α_Ratio", "IFN-γ_Ratio"),
    names_to = "Cytokine",
    values_to = "Ratio"
  )
df_summary <- df_long %>%
  group_by(Cell_Group, Cytokine) %>%
  summarise(
    mean = mean(Ratio),
    sd = sd(Ratio),
    .groups = "drop"
  )

library(cowplot)
plot_fun <- function(cyt){
  
  df_plot <- df_long %>% filter(Cytokine == cyt)
  df_sum  <- df_summary %>% filter(Cytokine == cyt)
  
  # 自定义颜色
  my_colors <- c(
    "#1f77b4", # 蓝
    "#ff7f0e", # 橙
    "#2ca02c", # 绿
    "#d62728", # 红
    "#9467bd", # 紫
    "#8c564b", # 棕
    "#e377c2", # 粉
    "#7f7f7f", # 灰
    "#bcbd22", # 黄绿
    "#17becf"  # 青
  )
  
  ggplot() +
    # bar (mean)
    geom_col(
      data = df_sum,
      aes(x = Cell_Group, y = mean),
      fill = "white",
      color = "black",
      width = 0.7
    ) +
    # SD error bar
    geom_errorbar(
      data = df_sum,
      aes(x = Cell_Group, ymin = mean - sd, ymax = mean + sd),
      width = 0.2,
      size = 0.6
    ) +
    # points (individual samples), now color = Cell_group
    geom_point(
      data = df_plot,
      aes(x = Cell_Group, y = Ratio, color = Cell_Group),
      fill = "white",
      size = 2.5,
      shape = 21,
      stroke = 1,
      position = position_jitter(width = 0.15)
    ) +
    scale_color_manual(values = my_colors) +
    theme_cowplot() +   # 使用 cowplot 主题
    labs(
      x = "",
      y = cyt,
      # title = paste("Expression:", cyt),
      color = "Cell group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

p_IL2  <- plot_fun("IL-2_Ratio")
p_TNF  <- plot_fun("TNF-α_Ratio")
p_IFNG <- plot_fun("IFN-γ_Ratio")

p_IL2
p_TNF
p_IFNG

stat_res <- df_long %>%
  group_by(Cytokine) %>%       # 每个细胞因子分别检验
  wilcox_test(Ratio ~ Cell_Group) %>% 
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat_res

(p_IL2 | p_TNF | p_IFNG) + plot_layout(guides = "collect")

# fig3h
# mki67+的CD8+T细胞大部分是tex-eff细胞
FeaturePlot(cd8t_tmp, "MKI67")
DimPlot(cd8t_tmp, group.by = "seurat_clusters") + ggsci::scale_color_igv()
DimPlot(cd8t_tmp, group.by = "cell_class") + ggsci::scale_color_igv()
cd8t_tmp$free_annotation <- ifelse(cd8t_tmp$seurat_clusters %in% c(9,10), "MKI67pos", "MKI67neg")
tmp_info <- dplyr::select(cd8t_tmp@meta.data, "free_annotation", "cell_class")
tmp_info$cell_class <- ifelse(tmp_info$cell_class == "Tex_eff", "Tex_eff", "else")
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = c("else","Tex_eff"))
ratio <- tmp_info %>%
  dplyr::count(free_annotation, cell_class) %>% 
  group_by(free_annotation) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n))

ggplot(data = ratio, aes(x=free_annotation, y=proportion, fill=cell_class)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="", fill="Cell Class") +
  coord_flip() +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = c("Tex_eff" = "#f9d580",
                               "else" = "lightgrey"))

# fig3i
tmp_info <- dplyr::select(cd8t_tmp@meta.data, "cell_class", "Classification","exclusion", "CTaa")
tmp_info <- na.omit(tmp_info) %>% dplyr::filter(!exclusion == "excluded")
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = c("Tnaive","Ttsm","Tpex","Tex_int","Tex_term","Tex_eff"))
tmp_info %>%
  dplyr::count(cell_class, Classification) %>% 
  group_by(cell_class) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=cell_class, y=proportion, fill=Classification)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Cell Class", fill="Propotion") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = c("else" = "#dddce7",
                               "Tumor Specific" = "#c571ac"))

# fig3j
# cd8t_tmp，Tex-eff中大克隆的占比最多
tmp_info <- dplyr::select(cd8t_tmp@meta.data, "cell_class", "cloneType", "sampleType")
tmp_info <- na.omit(tmp_info)
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = c("Tnaive","Ttsm","Tpex","Tex_int","Tex_term","Tex_eff"))
tmp_info %>%
  dplyr::count(cell_class, cloneType) %>% 
  group_by(cell_class) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=cell_class, y=proportion, fill=cloneType)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Cell Class", fill="Propotion") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  ggsci::scale_fill_igv()

in.dat <- dplyr::select(cd8t_tmp@meta.data, CTaa, cell_class, sampleType, PID, clone.status)
in.dat$Cell_Name <- rownames(cd8t_tmp@meta.data)
head(in.dat)
colnames(in.dat) <- c("clone.id", "majorCluster", "loc", "patient", "clone.status", "Cell_Name")
tic("Startrac.run")
out <- Startrac.run(in.dat, verbose=F)
toc()
Startrac::plot(out,index.type="cluster.all", byPatient=F)

# fig3k
head(cd8t_tmp@meta.data)
tmp_info <- dplyr::select(cd8t_tmp@meta.data, "NeoTCR8", "cell_class")
compaired <- list(c("Tnaive","Ttsm"),
                  c("Ttsm","Tpex"),
                  c("Tpex","Tex_int"),
                  c("Tex_term","Tex_int"),
                  c("Tex_eff","Tex_term"))
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = c("Tnaive","Ttsm","Tpex","Tex_int","Tex_term","Tex_eff"))
ggplot(data = tmp_info, aes(x = cell_class, y = NeoTCR8)) +
  geom_violin(aes(fill = cell_class), scale = "width", color = "white") +
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +
  stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +
  scale_x_discrete(limits = c("Tnaive", "Ttsm", "Tpex", "Tex_int", "Tex_term", "Tex_eff")) +  # 手动设置顺序
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = ex_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig3l
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = free_annotation)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("CD8_Texeff" = "#f9d580",
                                "Tumor reactive cell" = "#57afc6")) + 
  theme_classic() +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

# fig3m
cd8t_1$free_annotation <- ifelse(cd8t_1$cell_annotation == "CD8T_CXCL13", "CXCL13+CD8+T", "CXCL13-CD8+T")
cd8t_1 <- SetIdent(cd8t_1, value = "free_annotation")
# 提取UMAP坐标和注释数据
umap_data <- FetchData(cd8t_1, vars = c("UMAP_1", "UMAP_2", "free_annotation"))
# CD8T 细胞配色
cd8t_colors <- c("CXCL13+CD8+T" = "#166a3b","CXCL13-CD8+T" = "lightgrey")
# 绘制UMAP图
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = free_annotation)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = cd8t_colors) + 
  theme_classic() +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

# fig3n
tmp_info <- dplyr::select(cd8t_1@meta.data, "cell_class", "free_annotation", "sampleType")
tmp_info <- na.omit(tmp_info)
tmp_info$cell_annotation <- factor(tmp_info$cell_annotation, 
                                   levels = c("CXCL13+CD8+T",
                                              "CXCL13-CD8+T"))
tmp_info$cell_class <- factor(tmp_info$cell_class, c("Tnaive", "Ttsm", "Tpex", "Tex_int", "Tex_term", "Tex_eff"))
tmp_info %>%
  dplyr::count(cell_class, free_annotation) %>% 
  group_by(cell_class) %>%           
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=cell_class, y=proportion, fill=free_annotation)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Cell State", fill="Proportion") +
  theme_classic() +
  scale_fill_manual(values = cd8t_colors) +
  theme(legend.position = "none")  # 隐藏 legend

# fig4a
ggplot(data = cd8t_tmp.df[cd8t_tmp.df$cell_class %in% c("Tnaive","Ttsm","Tpex") & cd8t_tmp.df$State = "MSI" & cd8t_tmp.df$sampleType = "PLN",], aes(x = naive.score, y = prog.score)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.8, bins = 7, contour = TRUE) +  # 使用after_stat修正，调整透明度
  scale_fill_gradientn(colors = c("#fcf1f0", "#fccccb", "#e3716e"),  # 设置颜色渐变
                       values = scales::rescale(c(0, 0.5, 1)),  # 调整颜色在密度级别上的位置
                       limits = c(0, NA)) +  # 设置下限，使低密度区域有颜色
  geom_hline(yintercept = 0.14, linetype = "dashed", color = "black", size = 0.5) +
  geom_segment(aes(y = 0.14, yend = max(cd8t_tmp.df$prog.score), x = 0.3, xend = 0.3), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  scale_y_continuous(name = "Memory score") +
  scale_x_continuous(name = "Naive score")

ggplot(data = cd8t_tmp.df[cd8t_tmp.df$cell_class %in% c("Tnaive","Ttsm","Tpex") & cd8t_tmp.df$State = "MSS" & cd8t_tmp.df$sampleType = "PLN",], aes(x = naive.score, y = prog.score)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.8, bins = 7, contour = TRUE) +  # 使用after_stat修正，调整透明度
  scale_fill_gradientn(colors = c("#fcf1f0", "#fccccb", "#e3716e"),  # 设置颜色渐变
                       values = scales::rescale(c(0, 0.5, 1)),  # 调整颜色在密度级别上的位置
                       limits = c(0, NA)) +  # 设置下限，使低密度区域有颜色
  geom_hline(yintercept = 0.14, linetype = "dashed", color = "black", size = 0.5) +
  geom_segment(aes(y = 0.14, yend = max(cd8t_tmp.df$prog.score), x = 0.3, xend = 0.3), linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() +
  scale_y_continuous(name = "Memory score") +
  scale_x_continuous(name = "Naive score")

# fig4b
LN_data <- LN_data %>%
  mutate(X = (XMin + XMax) / 2) %>%
  mutate(Y = (YMin + YMax) / 2)
LN_data_info <- dplyr::select(LN_data, "Analysis Region", "CD8+PD1+TCF1+", "X", "Y")
LN_data_info$cell_annotation <- ifelse(LN_data_info$`CD8+PD1+TCF1+` == 1, "CD8_Tpex", "else")
color_map <- c(
  "else" = "#EAEAEA",     
  "CD8_Tpex" = "#bdb5e1"
)

ggplot() +
  # 先画“else”层
  geom_point(
    data = subset(LN_data_info, cell_annotation == "else"),
    aes(x = X, y = Y, color = cell_annotation),
    size = .1
  ) +
  geom_point(
    data = subset(LN_data_info, cell_annotation %in% c("CD8_Tpex")),
    aes(x = X, y = Y, color = cell_annotation),
    size = .2
  ) +
  # 自定义颜色映射
  scale_color_manual(values = color_map) +
  # 自定义图例点大小
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_cowplot() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# fig4c
tmp_info <- data.frame(cd8t_tmp[, cd8t_tmp$sampleType== "CRT"]$PID, cd8t_tmp[, cd8t_tmp$sampleType == "CRT"]$cell_class)
names(tmp_info) <- c("PID","seurat_clusters")
tmp_info <- group_by(tmp_info, PID, seurat_clusters) %>%
  summarise(n = n())
aa <- group_by(tmp_info, PID) %>%
  summarise(SUM = sum(n))

tmp_info <- merge(tmp_info, aa, by = "PID")
tmp_info$pct <- tmp_info$n/tmp_info$SUM
tmp_info$State <- ifelse(tmp_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")

ms_colors <- c("MSI" = "#ED0000FF",
               "MSS" = "#00468BFF")
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info[tmp_info$seurat_clusters %in% c("Tex_int","Tex_eff","Tex_term"),], 
       aes(x = State, y = pct, color = State)) +
  geom_text(aes(label = PID)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  scale_y_continuous(name = "Proportion of CD8+T cells in CRT") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  facet_wrap(~seurat_clusters, scales = "free_y") +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig4e
T_data <- T_data %>%
  mutate(X = (XMin + XMax) / 2) %>%
  mutate(Y = (YMin + YMax) / 2)
T_data_info <- dplyr::select(T_data, "Analysis Region", "CD8+GZMB+TIM3+",
                                    "X", "Y")
T_data_info$cell_annotation <- ifelse(T_data_info$`CD8+GZMB+TIM3+` == 1, "CD8_Texeff", T_data_info$`Analysis Region`)
table(T_data_info$cell_annotation)
color_map <- c(
  "IM" = "#acd6ea",       
  "NR" = "#cfe0a6",
  "TC" = "#fcf1f0",
  "CD8_Texeff" = "#f9d580"
)

ggplot() +
  geom_point(
    data = subset(T_data_info, !cell_annotation == "CD8_Texeff"),
    aes(x = X, y = Y, color = cell_annotation),
    size = .1
  ) +
  geom_point(
    data = subset(T_data_info, cell_annotation %in% c("CD8_Texeff")),
    aes(x = X, y = Y, color = cell_annotation),
    size = .2
  ) +
  # 自定义颜色映射
  scale_color_manual(values = color_map) +
  # 自定义图例点大小
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_cowplot() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# fig4f
compaired <- c("MSI","MSS")
ggplot(data = cd8texeff_info, aes(x = State, y = Ratio_CD8, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Ratio of CD8_Texeff/CD8+T") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  facet_wrap(.~Sample) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

cd8texeff_info$State <- factor(cd8texeff_info$State)
oneway_test(Ratio_CD8 ~ State, data = cd8texeff_info[cd8texeff_info$Sample == "Normal",], distribution = "exact") 

# fig4g
# 定义颜色
colors_late <- c("#FDAF91FF", adjustcolor("#FDAF91FF", alpha.f = 0.2))
colors_early <- c("#925E9FFF", adjustcolor("#925E9FFF", alpha.f = 0.2))

# 提取每组的比例
late_values <- table_data[table_data$Stage == "late", c("high", "low")]
early_values <- table_data[table_data$Stage == "early", c("high", "low")]

# 绘制饼图
par(mfrow = c(1, 2))  # 创建1行2列的布局

# Late组饼图
pie(
  as.numeric(late_values),
  labels = c("High", "Low"),
  col = colors_late,
  main = "Late Group",
  border = "white"  # 设置边界颜色为白色
)

# Early组饼图
pie(
  as.numeric(early_values),
  labels = c("High", "Low"),
  col = colors_early,
  main = "Early Group",
  border = "white"  # 设置边界颜色为白色
)

# 恢复默认布局
par(mfrow = c(1, 1))

compaired <- list(c("early","late"))
ggplot(data = tmp_info, aes(x = Stage, y = pct, color = Stage)) +
  geom_boxplot(aes(fill = Stage), alpha = 0.2, outlier.shape = NA) +  # 避免与抖动点重叠
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +  # 添加横向抖动
  scale_y_continuous(name = "Pct of CD8_Texeff cells") +
  scale_x_discrete(name = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.position = "none"
  ) +
  scale_color_manual(values = stage_color) +
  scale_fill_manual(values = stage_color)

# fig4h
aa <- tmp_info[tmp_info$Stage == "early",]
ggplot(data = aa, aes(x = State, y = pct, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Pct of CD8_Texeff cells") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) 

aa <- tmp_info[tmp_info$Stage == "late",]
ggplot(data = aa, aes(x = State, y = pct, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Pct of CD8_Texeff cells") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors)

tmp_info$State <- factor(tmp_info$State)
aa <- tmp_info[tmp_info$Stage == "early",] 
aa <- tmp_info[tmp_info$Stage == "late",] 
oneway_test(pct ~ State, data = aa, distribution = "exact") 

# fig4i
library(GSVA)
score_list <- list("tex.eff.score" = tex_eff.gene)
coad_matrix.tmp <- as.matrix(coad_matrix)
ssgsea_coad_matrix <- gsva(coad_matrix.tmp, score_list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
ssgsea_coad_matrix <- as.data.frame(t(ssgsea_coad_matrix))

ssgsea_coad_matrix <- ssgsea_coad_matrix[!duplicated(str_sub(rownames(ssgsea_coad_matrix),1,12)),,drop = FALSE]
rownames(ssgsea_coad_matrix) <- str_sub(rownames(ssgsea_coad_matrix),1,12)
ssgsea_coad_matrix$barcode <- rownames(ssgsea_coad_matrix) #453

head(ssgsea_coad_matrix)
head(survival_data)
head(clinical.coad.tmp)
survival_data <- arrange(survival_data, desc(futime)) #573
survival_data_tmp <- survival_data[!duplicated(survival_data$barcode),] #421

tex.eff_info <- merge(ssgsea_coad_matrix, clinical.coad.tmp, by.x = "barcode", by.y = "bcr_patient_barcode")
tex.eff_info.tmp <- dplyr::filter(tex.eff_info, !is.na(stage))

Stage_color <- c("early stage" = "#925E9FFF",
                 "late stage" = "#FDAF91FF")
tex.eff_info$Stage <- ifelse(tex.eff_info$stage %in% c("Stage I","Stage II"), "early stage",
                             ifelse(tex.eff_info$stage %in% c("Stage III", "Stage IV"), "late stage", NA))
tex.eff_info.tmp <- dplyr::filter(tex.eff_info, !is.na(stage))
head(tex.eff_info.tmp)
tex.eff_info.tmp <- merge(tex.eff_info.tmp, msi_info[,c(1,3)], by.x = "barcode", by.y = "PID")

compaired <- list(c("MSS","MSI"))
ggplot(data = tex.eff_info.tmp, aes(x = Status, y = tex.eff.score, color = Stage)) +
  geom_boxplot(aes(fill = Stage), alpha = 0.5, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1) +  # 添加横向抖动
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "CD8.Texeff score") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11)) +
  scale_color_manual(values = Stage_color) +
  scale_fill_manual(values = Stage_color) +  # 设置 fill 的颜色
  theme(legend.position = "none") +
  geom_signif(comparisons = compaired,
              step_increase = 0.3,
              map_signif_level = F,
              test = wilcox.test,
              color = "black")

# fig4j
library(survival)
library(survminer)
res.cut <- surv_cutpoint(tex.eff_info.tmp[tex.eff_info.tmp$Status == "MSS", ],
                         time = "futime", event = "fustat",
                         variables = "tex.eff.score")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(futime, fustat) ~ tex.eff.score, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = FALSE, pval = TRUE, 
           palette = c("#3c537a", "#ee9798"))

# fig4k
ggplot(ssgsea_exprSet_matrix.tmp, aes(x = Treat, y = tex.eff.score)) +
  geom_boxplot(aes(fill = Treat, color = Treat), alpha = 0.5) +
  geom_point(aes(color = Treat, group = PID)) +
  geom_line(aes(group = PID), color = "grey", linetype = "dotted") +
  theme_classic() +
  scale_fill_manual(values = c("#3c537a", "#ee9798")) +  # 设置箱线图填充颜色
  scale_color_manual(values = c("#3c537a", "#ee9798")) +  # 设置点和线的颜色
  labs(x = "Treatment", y = "Tex_eff Score") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11))

# fig4l
compaired <- list(c("response","non-response"))
ggplot(data = ssgsea_exprSet_matrix[ssgsea_exprSet_matrix$Treat == "Pre_xRT",], aes(x = Class, y = tex.eff.score, color = Class)) +
  geom_boxplot(aes(fill = Class), alpha = 0.5) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "CD8.Texeff score") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11)) +
  scale_color_manual(values = response_color) +
  scale_fill_manual(values = response_color) +  # 设置 fill 的颜色
  # facet_wrap(~seurat_clusters) +
  geom_signif(comparisons = compaired,
              step_increase = 0.3,
              map_signif_level = F,
              test = wilcox.test,
              color = "black")

# fig4m
compaired <- list(c("CR", "PR"),
                  c("CR", "SD"),
                  c("PR", "SD"))
response_color <- c("CR" = "#8C9EBA",
                    "PR" = "#D9E0E7", 
                    "SD" = "#F3DBD6")
ggplot(data = df[df$treat == "Pre",], aes(x = Response, y = pct, fill = Response)) +
  geom_boxplot() +
  # geom_point() +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Pct of CD8T in CRT")+
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  scale_fill_manual(values = response_color) +
  facet_wrap(~seurat_clusters, scales = "free_y") +
  geom_signif(comparisons = compaired,
              step_increase = 0.1,
              map_signif_level = F,
              test = wilcox.test)

# fig4n
library(survival)
library(survminer)
res.cut <- surv_cutpoint(tex.eff_info.tmp,
                         time = "futime", event = "fustat",
                         variables = "simple.score")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(futime, fustat) ~ simple.score, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = FALSE, pval = TRUE, 
           palette = c("#3c537a", "#ee9798"))

# fig4o
compaired <- list(c("response","non-response"))
ggplot(data = ssgsea_dat.tmp, aes(x = respond, y = simple.score, color = respond)) +
  geom_boxplot(aes(fill = respond), alpha = 0.5) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "IP score") +
  scale_x_discrete(name = "") +
  # ggtitle("Boxplot of CD8+ Tcell") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11)) +
  scale_color_manual(values = response_color) +
  scale_fill_manual(values = response_color) +  # 设置 fill 的颜色
  # facet_wrap(~seurat_clusters) +
  geom_signif(comparisons = compaired,
              step_increase = 0.3,
              map_signif_level = F,
              test = wilcox.test,
              color = "black")

# fig4p
immdat.exp_matrix.1 <- gsva(immdat.exp, score_list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
immdat.exp_matrix.1 <- as.data.frame(t(immdat.exp_matrix.1))
immdat.exp_matrix.1 <- merge(immdat.exp_matrix.1, p_inf[,c("Cancer","pfse","pfs","respond","best.resp")], by = 0)
str(immdat.exp_matrix.1)
immdat.exp_matrix.1$pfs.time <- as.numeric(immdat.exp_matrix.1$pfs)
immdat.exp_matrix.1$pfs.event <- as.numeric(immdat.exp_matrix.1$pfse)
immdat.exp_matrix.1$respond.event <- ifelse(immdat.exp_matrix.1$respond == "respond", 1, 0)

#建立曲线
rocobj1 <- roc(immdat.exp_matrix.1$respond.event, immdat.exp_matrix.1$tex.eff.score)
rocobj2 <- roc(immdat.exp_matrix.1$respond.event, immdat.exp_matrix.1$simple.score)
rocobj3 <- roc(immdat.exp_matrix.1$respond.event, immdat.exp_matrix.1$mhc.score)
#计算full AUC
auc(rocobj1)
auc(rocobj2)
auc(rocobj3)

# 绘制第一个 ROC 曲线
plot.roc(rocobj1, col = "black", lwd = 2)
# 添加第二个 ROC 曲线
lines(rocobj2, col = "darkred", lwd = 2)
# 添加第三个 ROC 曲线
lines(rocobj3, col = "grey", lwd = 2)
# 添加图例
legend("bottomright", 
       legend = c(paste("CD8.Texeff score (AUC =", round(auc(rocobj1), 2), ")"),
                  paste("IP score (AUC =", round(auc(rocobj2), 2), ")"),
                  paste("CD74 (AUC =", round(auc(rocobj3), 2), ")")),
       col = c("black", "darkred","grey"), 
       lwd = 1, 
       bty = "n",     # 去掉图例框
       cex = 0.8)     # 缩小图例大小

# fig4q
library(survival)
library(survminer)
res.cut <- surv_cutpoint(immdat.exp_matrix,
                         time = "pfs.time", event = "pfs.event",
                         variables = "simple.score",
                         minprop = 0.2)
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(pfs.time, pfs.event) ~ simple.score, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = FALSE, pval = TRUE, 
           palette = c("#3c537a", "#ee9798"))

# fig4r
response_color <- c("response" = "#9ab1d6",
                    "non-response" = "#fbc690")
immdat.exp_matrix$respond <- ifelse(immdat.exp_matrix$respond.event == 1, "response", "non-response")
ggplot(immdat.exp_matrix, aes(x = class, fill = respond)) +
  geom_bar(position = "fill") +  # 堆积条形图，显示比例
  scale_y_continuous(labels = scales::percent) +  # Y轴显示为百分比
  labs(x = "", y = "Percentage", fill = "Response",
       title = "") +
  scale_fill_manual(values = response_color) +
  theme_classic()

# fig4s
res.cut <- surv_cutpoint(dat,
                         time = "PFS", event = "PFS_event",
                         variables = "RiskScore")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(PFS, PFS_event) ~ RiskScore, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = FALSE, pval = TRUE, 
           palette = c("#3c537a", "#ee9798"))

# fig4t
t1 <- 1825
roc_3gene <- timeROC(
  T = dat$PFS,
  delta = dat$PFS_event,
  marker = dat$RiskScore,
  cause = 1,
  times = t1,
  iid = TRUE
)

roc_cxcl13 <- timeROC(
  T = dat$PFS,
  delta = dat$PFS_event,
  marker = dat$RiskScore_cxcl13,
  cause = 1,
  times = t1,
  iid = TRUE
)

auc_3gene  <- roc_3gene$AUC
auc_cxcl13 <- roc_cxcl13$AUC

# 5) 画在同一张图：先画三基因，再 add=TRUE 叠加 CXCL13
plot(roc_3gene, time = t1, col = "#1F77B4", lwd = 2, title = FALSE)
plot(roc_cxcl13, time = t1, col = "#D62728", lwd = 2, add = TRUE)

legend(
  "bottomright",
  legend = c("IP score", "CXCL13"),
  col    = c("#1F77B4","#D62728"),
  lwd    = 2,
  bty    = "n"
)

# fig5a
CD4T_colors <- c(
  "CD4T_CCR7" = "#a1c7d9",
  "CD4T_ANXA2_Tm" = "#f5b364",
  "CD4T_FOXP3" = "#eb8792",
  "CD4T_GZMK" = "#c18a6c",
  "CD4T_CXCL13" = "#684797",
  "CD4T_CXCL13/CXCR5" = "#3d659e",
  "CD4T_CXCR5" = "#7ba5be",
  "CD4T_GZMA" = "#e41e25",
  "CD4T_IL17A" = "#BDA7CB",
  "CD4T_ISG15" = "#166a3b",
  "CD4T_ANXA1/BHLHE40" = "#66a855",
  "CD4T_ANXA1_Tm" = "#add375"
)

DimPlot(cd4t_2, reduction = "umap", group.by = "cell_annotation", label = F, raster = FALSE, cols = color_palette) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)

# fig5b
cd4t_colors <- c("CD4T_Texeff" = "#8357ad",
                 "CD4T_Texterm" = "#759aae",
                 "CD4T_Texint" = "#b1cfea",
                 "CD4T_Tpex" = "#bbde92",
                 "CD4T_Tnaive" = "#8ebda7",
                 "Treg" = "#eb8792",
                 "Th17" = "#BDA7CB",
                 "Tfh" = "#3d659e")

DimPlot(cd4t_2, reduction = "umap", group.by = "cell_class", label = F, raster = FALSE, cols = cd4t_colors) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)

# fig5c
# 不同CD4亚群的特征基因表达
cd4t_tmp <- ScaleData(cd4t_tmp, features = rownames(cd4t_tmp))
Idents(cd4t_tmp) <- cd4t_tmp$cell_class
CD4T_markers <- c("SELL","LEF1","TCF7","CCR7","IL7R",
                  "PRF1","GZMB","GZMK","GZMA","NKG7","GNLY","IFNG","EOMES","TBX21",
                  "LAG3","TIGIT","PDCD1","CTLA4","ENTPD1","LAYN","HAVCR2","TOX",
                  "ICOS","CD28","CD40LG","TNFRSF14","TNFRSF9","TNFRSF18","TNFRSF4",
                  "CD200","CXCR5","CXCL13","BCL6","IL21",
                  "IL17A","IL17F","IL26","IL23R","IL4I1","RORC",
                  "IL2RA","FOXP3","CCR8","IKZF2","IKZF3",
                  "CD69","ITGAE","BHLHE40","ANXA1","ANXA2")
cd4t_matrix <- AverageExpression(cd4t_tmp, slot = "scale.data")[[1]]  
order_index <- match(CD4T_markers, rownames(cd4t_matrix))
col_index <- match(c("CD4T_Tnaive","CD4T_Tpex","CD4T_Texint","CD4T_Texeff","CD4T_Texterm","Tfh","Th17","Treg"),colnames(cd4t_matrix))
cd4t_matrix <- cd4t_matrix[order_index, ]
cd4t_matrix <- cd4t_matrix[, col_index]
cd4t_matrix[cd4t_matrix > 2] = 2;cd4t_matrix[cd4t_matrix < -2] = -2

pheatmap::pheatmap(cd4t_matrix, scale = "row", cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7", "#F3DBD6", "#DA8F87", "#D54846"))(100),
                   gaps_row = c(5,14,22,29,34,40,45),
                   gaps_col = c(1,2,3,4,5,6,7),
                   border_color = "white"
)

# fig5d
tmp_info <- select(cd4t_tmp@meta.data, "sampleType", "cell_class")
tmp_info$sampleType <- factor(tmp_info$sampleType, levels = c("PLN", "WBC", "CRN", "CRT"))
tmp_info$cell_class <- factor(tmp_info$cell_class, levels = rev(c("CD4T_Tnaive","CD4T_Tpex","CD4T_Texint","CD4T_Texeff","CD4T_Texterm","Tfh","Th17","Treg")))
tmp_info %>%
  dplyr::count(sampleType, cell_class) %>% 
  group_by(sampleType) %>%           # 按sampleType分组
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x=sampleType, y=proportion, fill=cell_class)) +  # 创建堆积条形图
  geom_bar(stat="identity") +
  labs(y="Proportion", x="Sample Type", fill="Cell Class") +
  theme_classic() +
  # ggtitle("Proportional Stacked Bar Chart of Cell Class by Sample Type") +
  scale_fill_manual(values = cd4t_colors)

# fig5f-g
tmp_info <- data.frame(cd4t_tmp$PID, cd4t_tmp$orig.ident, cd4t_tmp$cell_class)
names(tmp_info) <- c("PID","SID","seurat_clusters")
tmp_info <- dplyr::group_by(tmp_info, SID, seurat_clusters) %>%
  dplyr::summarise(n = n())
aa <- dplyr::group_by(tmp_info, SID) %>%
  dplyr::summarise(SUM = sum(n))
tmp_info <- merge(tmp_info, aa, by = "SID")
tmp_info$pct <- tmp_info$n/tmp_info$SUM

tmp_info$State <- ifelse(str_sub(tmp_info$SID,1,3) %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
tmp_info.1 <- filter(tmp_info, seurat_clusters %in% c("Treg","CD4T_Tpex"))
tmp_info.1$sampleType <- str_sub(tmp_info.1$SID,5,7)
tmp_info.1 <- dplyr::filter(tmp_info.1, sampleType == "PLN")
tmp_info.1$State <- factor(tmp_info.1$State, levels = c("MSI","MSS"))
tmp_info.1 <- arrange(tmp_info.1, desc(State))

ggplot(data = tmp_info.1,  aes(x = State, y = pct, color = State)) +
  geom_violin(aes(fill = State), scale = "width", alpha = 0.2) +
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  facet_wrap(~seurat_clusters, scales = "free_y") +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

result <- tmp_info.1 %>%
  filter(seurat_clusters %in% c("CD4T_Tpex", "Treg")) %>%  # 仅选择 CD4T_Tpex 和 Treg
  group_by(SID, seurat_clusters) %>%  # 按 SID 和 seurat_clusters 分组
  summarise(n_sum = sum(n, na.rm = TRUE), .groups = "drop") %>%  # 计算每个分组的数量
  pivot_wider(names_from = seurat_clusters, values_from = n_sum) %>%  # 将数据转为宽格式
  mutate(ratio = ifelse(!is.na(`CD4T_Tpex`) & !is.na(Treg), `CD4T_Tpex` / Treg, NA))  # 计算比值

print(result)
result$State <- ifelse(str_sub(result$SID,1,3) %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
result$PID <- str_sub(result$SID,1,3)
compaired <- list(c("MSI","MSS"))
ggplot(data = result, 
       aes(x = State, y = ratio, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID), color = "black", position = position_jitter(width = 0.1, height = 0)) +
  scale_y_continuous(name = "CD4_Tpex/Treg in TdLN") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  # facet_wrap(~seurat_clusters) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig5h
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info[tmp_info$sampleType == "CRT", ], 
       aes(x = State, y = pct, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  scale_y_continuous(name = "CD4T cells in CRT") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  facet_wrap(~seurat_clusters) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig5i
T_data_cd4 <- T_data_cd4 %>%
  mutate(X = (XMin + XMax) / 2) %>%
  mutate(Y = (YMin + YMax) / 2)
T_data_cd4_info <- dplyr::select(T_data_cd4, "Analysis Region", "CD4+IL21+IFNy+",
                                    "X", "Y")
T_data_cd4_info$cell_annotation <- ifelse(T_data_cd4_info$`CD4+IL21+IFNy+` == 1, "CD4_Texeff", T_data_cd4_info$`Analysis Region`)
table(T_data_cd4_info$cell_annotation)
color_map <- c(
  "IM" = "#acd6ea",       # 比 lightgrey 更浅的灰色
  "NR" = "#cfe0a6",
  "TC" = "#fcf1f0",
  "CD4_Texeff" = "#8357ad"
)

ggplot() +
  geom_point(
    data = subset(T_data_cd4_info, !cell_annotation == "CD4_Texeff"),
    aes(x = X, y = Y, color = cell_annotation),
    size = .1
  ) +
  geom_point(
    data = subset(T_data_cd4_info, cell_annotation %in% c("CD4_Texeff")),
    aes(x = X, y = Y, color = cell_annotation),
    size = .2
  ) +
  # 自定义颜色映射
  scale_color_manual(values = color_map) +
  # 自定义图例点大小
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_cowplot() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# fig5j-k
compaired <- c("MSI","MSS")
ggplot(data = cd4texeff_info, aes(x = State, y = Density, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Density of CD4_Texeff cells/ROI") +
  scale_x_discrete(name = "") +
  facet_wrap(.~Sample) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

ggplot(data = cd4texeff_info, aes(x = State, y = Ration_CD4, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.shape = NA) +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  # geom_text(aes(label = PID)) +
  scale_y_continuous(name = "Ratio of CD4_Texeff/CD4+T") +
  scale_x_discrete(name = "") +
  facet_wrap(.~Sample) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig5l
cd4t_info <- data.frame(cd4t_2[, cd4t_2$sampleType== "CRT"]$PID, cd4t_2[, cd4t_2$sampleType == "CRT"]$cell_annotation)
names(cd4t_info) <- c("PID","seurat_clusters")
cd4t_info <- group_by(cd4t_info, PID, seurat_clusters) %>%
  summarise(n = n())
aa <- group_by(cd4t_info, PID) %>%
  summarise(SUM = sum(n))
cd4t_info <- merge(cd4t_info, aa, by = "PID")
cd4t_info$pct <- cd4t_info$n/cd4t_info$SUM

tmp_info <- data.frame(cd8t_tmp[, cd8t_tmp$sampleType== "CRT"]$PID, cd8t_tmp[, cd8t_tmp$sampleType == "CRT"]$cell_class)
names(tmp_info) <- c("PID","seurat_clusters")
tmp_info <- group_by(tmp_info, PID, seurat_clusters) %>%
  summarise(n = n())
aa <- group_by(tmp_info, PID) %>%
  summarise(SUM = sum(n))
tmp_info <- merge(tmp_info, aa, by = "PID")
tmp_info$pct <- tmp_info$n/tmp_info$SUM
tmp_info$State <- ifelse(tmp_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
tmp_info <- merge(tmp_info, cd4t_info, by = "PID")
tmp_info_1 <- filter(tmp_info, seurat_clusters.x == "Tex_eff" & seurat_clusters.y == "CD4T_CXCL13")

library(ggpmisc)
model <- lm(pct.x ~ pct.y, data = tmp_info_1)
r_squared <- summary(model)$adj.r.squared
p_value <- summary(model)$coefficients[2, 4]
ggplot(tmp_info_1, aes(x = pct.y, y = pct.x, color = State)) +
  geom_point() +  # 添加点
  theme_classic() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_y_continuous(name = "Proportion of CD8_Texeff") +
  scale_x_continuous(name = "Proportion of CD4_Texeff") +
  # 显示R²和p值，不区分State
  annotate("text", x = 0.05, y = 0.95, label = sprintf("adj.R² = %.3f, p = %.3g", round(r_squared, digits = 2), round(p_value, digits = 2)),
           size = 5, hjust = 0, vjust = 1, parse = FALSE) +
  scale_color_manual(values = c("MSI" = "#ED0000FF", "MSS" = "#00468BFF"))

# fig5m
DotPlot(cd4t_tmp[, cd4t_tmp$sampleType == "CRT"], features = c("IFNG", "IL21"), group.by = "State") + 
  coord_flip() +
  labs(
    title = "",  # 设置图的标题
    x = "",     # 设置x轴标签
    y = ""      # 设置y轴标签
  ) +
  scale_color_gradient2(low = "#4575b4", high = "#d73027") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

DotPlot(cd8t_tmp[, cd8t_tmp$sampleType == "CRT"], features = c("IL21R"), group.by = "State") + 
  coord_flip() +
  labs(
    title = "",  # 设置图的标题
    x = "",     # 设置x轴标签
    y = ""      # 设置y轴标签
  ) +
  scale_color_gradient2(low = "#4575b4", high = "#d73027") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# fig5n
library(ggplot2)
plot <- FeaturePlot(cd8t_tmp, features = "IL21R")
data <- plot$data
colnames(data) <- c("x", "y", "ident", "gene")

# 重新排序数据框，基于 gene 表达值，从小到大排序
data <- data[order(data$gene, decreasing = FALSE), ]

# 绘制图形
ggplot(data, aes(x = x, y = y)) +
  geom_point(size = .1,
             aes(color = gene), show.legend = TRUE) +
  scale_color_gradientn(
    colors = c("lightgrey", "#b3a3c8", "#907bb3", "#5e529a", "#304693"),
    na.value = "lightgrey", 
    limits = c(NA, max(data$gene, na.rm = TRUE))
  ) +
  geom_density_2d(data = data[!is.na(data$gene) & data$gene > median(data$gene), ],
                  aes(x = x, y = y), bins = 3, color = "black") +
  theme_bw() +
  ggtitle("") +
  labs(x = "", y = "", color = "IL21R") +  # 修改图例名称为 IL2R
  theme(
    panel.grid = element_blank(),            # 去掉背景横线和竖线
    axis.title = element_blank(),            # 去掉 X 和 Y 轴的标题
    axis.text = element_blank(),             # 去掉 X 和 Y 轴的文字
    axis.ticks = element_blank()             # 去掉 X 和 Y 轴的刻度
  )

plot <- FeaturePlot(cd8t_tmp, features = "IFNGR1")
data <- plot$data
colnames(data) <- c("x", "y", "ident", "gene")

# 重新排序数据框，基于 gene 表达值，从小到大排序
data <- data[order(data$gene, decreasing = FALSE), ]

# 绘制图形
ggplot(data, aes(x = x, y = y)) +
  geom_point(size = .1,
             aes(color = gene), show.legend = TRUE) +
  scale_color_gradientn(
    colors = c("lightgrey", "#b3a3c8", "#907bb3", "#5e529a", "#304693"),
    na.value = "lightgrey", 
    limits = c(NA, max(data$gene, na.rm = TRUE))
  ) +
  geom_density_2d(data = data[!is.na(data$gene) & data$gene > median(data$gene), ],
                  aes(x = x, y = y), bins = 3, color = "black") +
  theme_bw() +
  ggtitle("") +
  labs(x = "", y = "", color = "IFNGR1") +  # 修改图例名称为 IL2R
  theme(
    panel.grid = element_blank(),            # 去掉背景横线和竖线
    axis.title = element_blank(),            # 去掉 X 和 Y 轴的标题
    axis.text = element_blank(),             # 去掉 X 和 Y 轴的文字
    axis.ticks = element_blank()             # 去掉 X 和 Y 轴的刻度
  )

# fig5o
import_para_summary$mean_importance[is.na(import_para_summary$mean_importance)] <- 0
top10_df <- import_para_summary %>%
  filter(Predictor == target_of_interest) %>%        
  group_by(Target) %>%
  summarise(
    mean_importance = mean(mean_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance)) %>%               # 从大到小排序
  slice_head(n = 10) %>%                           # 取前 10 个
  mutate(
    Rank = row_number(),                           # 1~10
    Target = factor(Target, levels = Target)  # 按重要性顺序固定颜色 legend
  )
top10_df

top10_df <- top10_df %>%
  mutate(
    Target = factor(Target, levels = Target)  # 按重要性顺序
  )

p <- ggplot(
  top10_df,
  aes(x = Target, y = mean_importance, color = Target, size = mean_importance)
) +
  geom_point(alpha = 1) +
  scale_size_continuous(range = c(3, 10)) +
  theme_cowplot() +
  labs(
    x = "",
    y = "Mean importance",
    title = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_manual(values = c("Mac_NTM_THBS1" = "#ffaf75",
                                "Mac_IL10" = "#ffaf75",
                                "Mac_MMP9" = "#d58f5c",
                                "gdT_TRDV1" = "#9869aa",
                                "NK_TRDV1_like" = "#cae4f3",
                                "NK_SELL" = "#74daf3",
                                "Bcell_Memory" = "grey",
                                cd8t_colors,
                                CD4T_colors))

p

# fig6a
dc.deg <- FindMarkers(subset(dc, sampleType == "CRT"),
                      ident.1 = "MSI", ident.2 = "MSS",
                      group.by = "State", min.pct = 0.1,  
                      logfc.threshold = 0.1)
deg <- dc.deg
# 火山图
library(ggrepel)
library(scales)
# 设置抖动
jitter_height <- 0.1
jitter <- position_jitter(height = jitter_height, seed = 1) ####Create a jitter object for reproducible jitter
# 添加logp与pct
correct_value = min(deg$p_val_adj[deg$p_val_adj > 0])
deg$logp = -1* log10(ifelse(deg$p_val_adj > 0, deg$p_val_adj, correct_value))
deg$pct <- ifelse(deg$avg_log2FC > 0, deg$pct.1, deg$pct.2)
deg$gene_change <- ifelse((deg$avg_log2FC > 0.2 & deg$p_val_adj < 0.05), "Up",
                          ifelse((deg$avg_log2FC < -0.2 & deg$p_val_adj < 0.05), "Down", "Else"))
table(deg$gene_change)
dc.deg.tmp <- filter(dc.deg, p_val_adj < 0.05)
# 设置需要展示的基因名称
dc.deg.tmp <- dc.deg.tmp %>% 
  dplyr::arrange(desc(avg_log2FC))
# n <- 30
genes <- c("RELB","STAT1","OAS1","IRF1","IRF7","IRF8","ISG15",
           "CLEC4A","CLEC10A","CD1C",
           "HLA-A","HLA-B",
           "IL1A","CXCL9","CXCL10","CCL5","CCL17","CCL19",
           "CLEC9A","XCR1","BATF3",
           "CCR7",
           "CLEC4C","IL3RA","TCL1A", rownames(dc.deg.tmp)[c(1:5,(nrow(dc.deg.tmp)-4):nrow(dc.deg.tmp))])
deg$gene_name <- ifelse(rownames(deg) %in% genes, rownames(deg), NA)
p <- deg %>% 
  dplyr::arrange(desc(avg_log2FC)) %>%
  ggplot(aes(x = avg_log2FC, y = logp, color = gene_change)) +
  geom_point(aes(size = pct), position = jitter, alpha = 1) +  # 使用带有纵向抖动的 jitter
  geom_text_repel(aes(x = avg_log2FC, y = logp, label = gene_name),
                  color = "black", position = jitter, size = 3, fontface = "bold", 
                  max.overlaps = 50, segment.color = "#808080") +  # 设置线的颜色为灰色
  theme_bw(base_size = 15) +
  ggtitle(paste0("DEGs of DC in MSI CRC TME")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +  # log2FC = -0.3 和 0.3
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p = 0.01
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        text = element_text(face = "bold")) +
  scale_color_manual(values=c("Down"="#00468BFF", "Else"="grey", "Up"="#ED0000FF"))
p

# fig6b
# 通路富集
library(clusterProfiler)
library(org.Hs.eg.db)
deg.sig.up <- subset(dc.deg.tmp, avg_log2FC > 0.5)
ego_BP <- enrichGO(gene          = row.names(deg.sig.up),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
ego_bp <- ego_bp[order(ego_bp$Count, decreasing = TRUE),]
ego_bp$Description <- factor(ego_bp$Description, levels = ego_bp$Description)
ego_bp <- separate(data = ego_bp, 
                   col = GeneRatio,
                   into = c("GR1", "GR2"), 
                   sep = "/")
ego_bp <- mutate(ego_bp, 
                 GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
# library(forcats)
ggplot(ego_bp[1:12,], aes(x = GeneRatio,y= fct_reorder(Description, GeneRatio)))+#将term顺序按照GeneRatio进行排序
  geom_point(aes(size = Count,fill = p.adjust),
             shape = 21,
             color = 'black')+
  #修改气泡图颜色
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+
  #标题修改
  labs(title='GO Enrichment',
       y='GO term',
       x='GeneRatio')+
  guides(fill=guide_colorbar(reverse = T,order=1))+
  theme_bw() +
  theme(axis.text = element_text(face = "bold")) 

# fig6c
# MSI DC的IFN响应评分
tmp_info <- dplyr::select(dc@meta.data, "PID", "sampleType", "IFN_MHCII_Score","IFNA_Response_Score","IFNG_Response_Score")
compaired <- list(c("MSS","MSI"))
tmp_info$State <- ifelse(tmp_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42"), "MSI", "MSS")
tmp_info <- filter(tmp_info, sampleType == "CRT")
tmp_info <- pivot_longer(tmp_info, cols = c(4,5), names_to = "IFN_Response_Pathway", values_to = "IFN_Score")
ggplot(data = tmp_info, aes(x = State, y = IFN_Score, fill = State)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  scale_y_continuous(name = "IFN Response Score")+
  scale_x_discrete(name = "State") +
  # stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    # axis.text.x = element_blank(),  # 隐藏横轴文字
    # axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ms_colors) +
  facet_wrap(~IFN_Response_Pathway) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig6d
deg_dc.tdln <- FindMarkers(subset(dc, sampleType == "PLN"),
                           ident.1 = "MSI", ident.2 = "MSS",
                           group.by = "State")

deg_dc.tdln.sig <- deg_dc.tdln %>% filter(p_val_adj < 0.05)
deg_dc.tdln.up <- rownames(filter(deg_dc.tdln.sig, avg_log2FC > 0.5))
deg_dc.tdln.down <- rownames(filter(deg_dc.tdln.sig, avg_log2FC < -0.5))

# fig6e
library(clusterProfiler)
library(org.Hs.eg.db)
ego_BP <- enrichGO(gene          = deg_dc.tdln.up,
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
ego_bp <- ego_bp[order(ego_bp$Count, decreasing = TRUE),]
ego_bp <- separate(data = ego_bp, 
                   col = GeneRatio,
                   into = c("GR1", "GR2"), 
                   sep = "/")
ego_bp <- mutate(ego_bp, 
                 GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
ggplot(ego_bp, aes(x = GeneRatio, y = Description)) +
  geom_segment(aes(x = 0, xend = GeneRatio, y = Description, 
                   yend = Description, color = -log10(p.adjust)), size = 1.5) + # 调整线段粗细
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_y_discrete(labels = ego_bp_top10$Description_wrapped) + # 使用换行后的标签
  scale_color_gradient(low = "#fdc5b8", high = "#d85152") +
  labs(
    x = "GeneRatio",
    y = NULL,
    # title = "Lollipop Plot for Pathway Enrichment",
    color = "-log10(P.adjust)",
    size = "Count"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10)
  )

# fig6f
tmp_info <- dplyr::select(dc@meta.data, "IFNA_Response_Score", "IFNG_Response_Score", "STING_Score", "State", "sampleType", "cell_annotation")
ggplot(data = tmp_info[tmp_info$sampleType == "PLN",], aes(x = State, y = IFNG_Response_Score, fill = State)) +
  geom_violin(scale = "width", color = "white") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.alpha = 0, color = "black") +
  theme_classic() +
  # theme(
  #   plot.title = element_text(size = 14, face = "bold"),
  #   text = element_text(size = 10),
  #   axis.title = element_text(face = "bold"),
  #   axis.title.x = element_blank(),  # 去掉 X 轴的标题
  #   axis.text.x = element_blank(),  # 隐藏横轴文字
  #   axis.ticks.x = element_blank(),  # 隐藏横轴刻度
  #   legend.position = "none"  # 移除图例
  # ) +
  scale_fill_manual(values = ms_colors) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = TRUE,
    test = wilcox.test
  )

# fig6g
mac.deg <- FindMarkers(subset(mac, sampleType == "CRT"),
                       ident.1 = "MSI", ident.2 = "MSS",
                       group.by = "State", min.pct = 0.1,  
                       logfc.threshold = 0.1)
deg <- mac.deg
# 添加logp与pct
correct_value = min(deg$p_val_adj[deg$p_val_adj > 0])
deg$logp = -1* log10(ifelse(deg$p_val_adj > 0, deg$p_val_adj, correct_value))
deg$pct <- ifelse(deg$avg_log2FC > 0, deg$pct.1, deg$pct.2)
deg$gene_change <- ifelse((deg$avg_log2FC > 0.2 & deg$p_val_adj < 0.05), "Up",
                          ifelse((deg$avg_log2FC < -0.2 & deg$p_val_adj < 0.05), "Down", "Else"))
table(deg$gene_change)
mac.deg.tmp <- filter(mac.deg, p_val_adj < 0.05)
genes <- c("TNFAIP3","TNF",
           "CD86","HLA-A","HLA-B",
           "CXCL9","CXCL10","CXCL11","CXCL16")
deg$gene_name <- ifelse(rownames(deg) %in% genes, rownames(deg), NA)
p <- deg %>% 
  dplyr::arrange(desc(avg_log2FC)) %>%
  ggplot(aes(x = avg_log2FC, y = logp, color = gene_change)) +
  geom_point(aes(size = pct), position = jitter, alpha = 1) +  # 使用带有纵向抖动的 jitter
  geom_text_repel(aes(x = avg_log2FC, y = logp, label = gene_name),
                  color = "black", position = jitter, size = 3, fontface = "bold", 
                  max.overlaps = 50, segment.color = "#808080") +  # 设置线的颜色为灰色
  theme_bw(base_size = 15) +
  ggtitle(paste0("DEGs of Macrophage in MSI CRC TME")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +  # log2FC = -0.3 和 0.3
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p = 0.01
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        text = element_text(face = "bold")) +
  scale_color_manual(values=c("Down"="#00468BFF", "Else"="grey", "Up"="#ED0000FF"))
p

# fig6h
library(clusterProfiler)
library(org.Hs.eg.db)
deg.sig.up <- subset(mac.deg.tmp, avg_log2FC > 0.5)
ego_BP <- enrichGO(gene          = row.names(deg.sig.up),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
ego_bp <- ego_bp[order(ego_bp$Count, decreasing = TRUE),]
ego_bp$Description <- factor(ego_bp$Description, levels = ego_bp$Description)
ego_bp <- separate(data = ego_bp, 
                   col = GeneRatio,
                   into = c("GR1", "GR2"), 
                   sep = "/")
ego_bp <- mutate(ego_bp, 
                 GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
# library(forcats)
ggplot(ego_bp[1:12,], aes(x = GeneRatio,y= fct_reorder(Description, GeneRatio)))+#将term顺序按照GeneRatio进行排序
  geom_point(aes(size = Count,fill = p.adjust),
             shape = 21,
             color = 'black')+
  #修改气泡图颜色
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+
  #标题修改
  labs(title='GO Enrichment',
       y='GO term',
       x='GeneRatio')+
  guides(fill=guide_colorbar(reverse = T,order=1))+
  theme_bw() +
  theme(axis.text = element_text(face = "bold")) 

# fig6i
DimPlot(mac, reduction = "umap", group.by = "cell_annotation", label = F, raster = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggsci::scale_color_igv() +
  labs(title = NULL)

# fig6j
tmp_info <- data.frame(mac[, mac$sampleType== "CRT"]$PID, mac[, mac$sampleType == "CRT"]$cell_annotation)
names(tmp_info) <- c("PID","seurat_clusters")
tmp_info <- group_by(tmp_info, PID, seurat_clusters) %>%
  summarise(n = n())
aa <- group_by(tmp_info, PID) %>%
  summarise(SUM = sum(n))
tmp_info <- merge(tmp_info, aa, by = "PID")
tmp_info$pct <- tmp_info$n/tmp_info$SUM
tmp_info$State <- ifelse(tmp_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info[tmp_info$seurat_clusters == "Mac_CXCL10",], 
       aes(x = State, y = pct, color = State)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  scale_y_continuous(name = "Mac_CXCL10 cells in CRT") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  facet_wrap(~seurat_clusters) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig6k-l
netVisual_aggregate(cellchat_obj, signaling = pathways.show, 
                    layout = "circle", color.use = NULL, 
                    sources.use = Mac_CXCL10, targets.use = NULL, 
                    idents.use = NULL)
netAnalysis_contribution(cellchat_obj, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat_obj, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[c(1,2,5),]
netVisual_individual(cellchat_obj, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_bubble(
  cellchat_obj,
  sources.use = "Mac_CXCL10",
  targets.use = c("CD8_Texeff", "CD4_Texeff"),
  signaling = "CXCL",
  remove.isolate = FALSE
)

# fig6m
ggplot(data = meta_mac,
       aes(x = level3, y = prop_log10, fill = level3, colour = level3)) +
  geom_sina(alpha = 1, size = 0.5) +
  geom_point(
    data = q90_df,
    aes(x = level3, y = q90),
    inherit.aes = FALSE,
    color = "grey20",
    size = 2
  ) +
  geom_line(
    data = q90_df,
    aes(x = level3, y = q90, group = 1),
    inherit.aes = FALSE,
    color = "grey50",
    linewidth = 0.5
  ) +
  theme_cowplot() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),   # 跟你原来一样处理
    axis.text.x  = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = tissue_colors) +
  scale_color_manual(values = tissue_colors) +
  labs(
    y = "log10(Mac_CXCL10 RTCD)",
    x = NULL,
    title = ""
  )

# fig6n
# ST
Idents(mac) <- mac$cell_annotation
mac_cxcl10.deg <- FindMarkers(mac,
                              ident.1 = "Mac_CXCL10",
                              min.pct = 0.25,
                              logfc.threshold = 0.1)
mac_cxcl10.deg.filter <- mac_cxcl10.deg %>% filter(p_val_adj < 0.05)
mac_cxcl10.deg.filter <- mac_cxcl10.deg.filter[order(mac_cxcl10.deg.filter$avg_log2FC, decreasing = T),]
mac_cxcl10.degs <- rownames(mac_cxcl10.deg.filter)[1:50]
DimPlot(mac)

markerList2 <- list(
  "Tex_eff.score" = tex_eff.gene,
  "Mac_CXCL10.score" = mac_cxcl10.degs,
  "tls.score" = tls.gene.12)
cells_ranking <- AUCell_buildRankings(st_object.p67@assays$Spatial@data)
cells_auc <- AUCell_calcAUC(markerList2, cells_ranking, aucMaxRank = 0.1 * nrow(cells_ranking))
aa <- as.data.frame(t(cells_auc@assays@data$AUC))
st_object.p67 <- AddMetaData(st_object.p67, metadata = aa)

SpatialPlot(st_object.p67,
            features = "Tex_eff.score",
            image.alpha = 0,
            pt.size.factor = 2,
            stroke = NA)

SpatialPlot(st_object.p67,
            features = "Mac_CXCL10.score",
            image.alpha = 0,
            pt.size.factor = 2,
            stroke = NA) 

FeaturePlot(st_object.p67, features = "Mac_CXCL10.score", reduction='spatial') + 
  scale_y_reverse() + 
  scale_colour_viridis(option="inferno") + theme_void()

# fig6o
import_para_summary$mean_importance[is.na(import_para_summary$mean_importance)] <- 0
top10_df <- import_para_summary %>%
  filter(Predictor == target_of_interest) %>%        
  group_by(Target) %>%
  summarise(
    mean_importance = mean(mean_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance)) %>%               # 从大到小排序
  slice_head(n = 10) %>%                           # 取前 10 个
  mutate(
    Rank = row_number(),                           # 1~10
    Target = factor(Target, levels = Target)  # 按重要性顺序固定颜色 legend
  )
top10_df

top10_df <- top10_df %>%
  mutate(
    Target = factor(Target, levels = Target)  # 按重要性顺序
  )

p <- ggplot(
  top10_df,
  aes(x = Target, y = mean_importance, color = Target, size = mean_importance)
) +
  geom_point(alpha = 1) +
  scale_size_continuous(range = c(3, 10)) +
  theme_cowplot() +
  labs(
    x = "",
    y = "Mean importance",
    title = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_manual(values = c("Mac_NTM_THBS1" = "#ffaf75",
                                "Mac_IL10" = "#ffaf75",
                                "Mac_MMP9" = "#d58f5c",
                                "gdT_TRDV1" = "#9869aa",
                                "NK_TRDV1_like" = "#cae4f3",
                                "NK_SELL" = "#74daf3",
                                "Bcell_Memory" = "grey",
                                cd8t_colors,
                                CD4T_colors))
p

# fig7a
tmp_info <- data.frame(bpls[, bpls$sampleType== "CRT"]$PID, bpls[, bpls$sampleType == "CRT"]$antibody_class)
names(tmp_info) <- c("PID","seurat_clusters")
tmp_info <- group_by(tmp_info, PID, seurat_clusters) %>%
  summarise(n = n())
aa <- group_by(tmp_info, PID) %>%
  summarise(SUM = sum(n))
tmp_info <- merge(tmp_info, aa, by = "PID")
tmp_info$pct <- tmp_info$n/tmp_info$SUM
tmp_info$State <- ifelse(tmp_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info[tmp_info$seurat_clusters %in% c("IGHG1"),], 
       aes(x = State, y = pct, color = State)) +
  # geom_text(aes(label = PID)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  scale_y_continuous(name = "Proportion of plasma cells in CRT") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  facet_wrap(~seurat_clusters, scales = "free_y") +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig7c
compaired <- list(c("MSI","MSS"))
tls_info$State <- ifelse(tls_info$PID %in% c("P02", "P03", "P08", "P17", "P24", "P27", "P42", "P44", "P48"), "MSI", "MSS")
ggplot(data = tls_info, 
       aes(x = State, y = ratio_1, color = State)) +
  # geom_text(aes(label = PID)) +
  geom_boxplot(aes(fill = State), alpha = 0.2, outlier.color = "white") +  # 使用 aes 设置 fill 参数
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1) +
  scale_y_continuous(name = "TLS/mm2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 x 轴的题目
    axis.text.x = element_text(size = 11),
    legend.position = "none"         # 去掉图例
  ) +
  scale_color_manual(values = ms_colors) +
  scale_fill_manual(values = ms_colors) +  # 设置 fill 的颜色
  # facet_wrap(~seurat_clusters, scales = "free_y") +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test,
    color = "black"
  )

# fig7d
library(corrplot)
corrplot(correlation_matrix, 
         method = "color",       # 使用颜色填充单元格
         type = "lower",         # 只显示矩阵的下三角部分
         order = "hclust",       # 使用层次聚类对相关矩阵进行排序
         tl.pos = "l",           # 仅在左侧显示行名，隐藏列名
         tl.col = "black",       # 行名颜色为黑色
         col = colorRampPalette(c("navyblue", "white", "darkred"))(50) # 颜色渐变从蓝色到白色再到红色
)

# fig7e
netVisual_individual(cellchat_obj, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# fig7f
ggplot(data = result_df, aes(x = Spearman_Correlation, y = Median_Ratio, colour = cell_annotation)) +
  geom_point(size = 6) +
  theme_classic() +
  geom_text(aes(label = cell_annotation), color = "black", 
            position = position_jitter(height = 0.2),
            size = 3) +  # 添加扰动
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = c(cd8t_colors, CD4T_colors,
                                "aDC" = "#ae1f63",
                                "Mac_CXCL10" = "#ce3d32",
                                "Mac_SPP1" = "#d58f5c")) +
  theme(legend.position = "none")  # 去掉图例

# fig7g
library(ggpmisc)
ggplot(data = ssgsea_coad_matrix, aes(x = tex.eff.score, y = tls12.score)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_poly_eq(use_label(c("adj.R2", "p.value.label")),
               formula = y ~ x,  parse = TRUE,
               size = 5, #公式字体大小
               label.x = 0.05,  #位置 ，0-1之间的比例
               label.y = 0.95)
ggplot(data = ssgsea_coad_matrix, aes(x = cd4t.cxcl13.score, y = tls12.score)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_poly_eq(use_label(c("adj.R2", "p.value.label")),
               formula = y ~ x,  parse = TRUE,
               size = 5, #公式字体大小
               label.x = 0.05,  #位置 ，0-1之间的比例
               label.y = 0.95)

# fig7h
SpatialPlot(st_object.p67,
            features = "tls.score",
            image.alpha = 0,
            pt.size.factor = 2,
            stroke = NA) 

# fig7k
library(survival)
library(survminer)
res.cut <- surv_cutpoint(tex.eff_info,
                         time = "futime", event = "fustat",
                         variables = "tls12.score")
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(futime, fustat) ~ tls12.score, data = res.cat)
ggsurvplot(fit, data = res.cat, risk.table = FALSE, pval = TRUE, 
           palette = c("#3c537a", "#ee9798"))

# fig8a
plot <- DimPlot(epi.nmf, reduction = "umap", group.by = "seurat_clusters", label = F, raster = FALSE) +
  ggsci::scale_color_lancet() +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL)

# fig8b
# 绘制点图
H_long$Gene <- factor(H_long$Gene, levels = combined_vector)
library(ggplot2)
library(viridis)
ggplot(H_long, aes(x = Gene, y = Gene_module, size = Weight, color = Weight)) +
  geom_point(alpha = 1) +  
  scale_size(range = c(1,6)) +  # 调整点大小范围
  scale_color_viridis_c(option = "D", direction = -1) +  # 使用 Viridis 色板
  theme_classic() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签以避免重叠
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    axis.title.y = element_blank()  
  ) +
  labs(
    x = "Gene",
    y = "Gene Module",
    size = "Weight",
    color = "Weight"
    # title = "Gene Contributions to Modules"
  )

# fig8c
df <- select(epi.nmf@meta.data, "PID","Stress_Response_Module","Antigen_Presentation_Module","Cell_Proliferation_Module",
             "IFNγ_Response_Module","Mitochondrial_Detoxification_Module","Tumor_Marker_Module","Intestinal_Epithelial_Marker_Module1",
             "Intestinal_Epithelial_Marker_Module2")
result <- aggregate(. ~ PID, data = df, FUN = mean, na.rm = TRUE)
rownames(result) <- result[,1]
result <- result[,-1]
colnames(result) <- paste0("P",1:8)
result <- result[c("P03","P17","P24","P42","P25","P38","P41"),c("Module2","Module4","Module7","Module8","Module1","Module3","Module5","Module6")]
result <- data.frame(t(result))
anno <- data.frame(State = c(rep("MSI", 4), rep("MSS", 3)),
                   row.names = colnames(result))

# 绘制热图
pheatmap::pheatmap(result, 
                   scale = "row", 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = anno,
                   annotation_colors = colors,
                   gaps_row = 4,
                   gaps_col = 4,
                   color = colorRampPalette(c("#547297", "#8C9EBA", "#D9E0E7","#F3DBD6", "#DA8F87", "#D54846"))(100),
                   annotation_legend = FALSE)

# fig8d
ggplot(data = ego_bp[1:10,], aes(x = Description, y = Count)) + 
  geom_bar(stat="identity", width=0.8, aes(fill = p.adjust)) + 
  coord_flip() + xlab("") + ylab("") + 
  theme_classic(base_size = 12) +
  scale_fill_gradient(low = "#AD002AFF", high = "#ADB6B6FF") +
  scale_x_discrete(labels=function(x) str_wrap(x, width = 50)) +
  theme(axis.text = element_text(face = "bold"),
        legend.position = "none")  # 加粗坐标轴文字

# fig8e
# 在MSI与MSS肿瘤上皮中找差异基因
epi.deg <- FindMarkers(epi.nmf,
                       ident.1 = "MSI",
                       ident.2 = "MSS",
                       group.by = "State")
deg <- epi.deg
# 设置抖动
jitter_height <- 0.1
jitter <- position_jitter(height = jitter_height, seed = 1) ####Create a jitter object for reproducible jitter
# 添加logp与pct
correct_value = min(deg$p_val_adj[deg$p_val_adj > 0])
deg$logp = -1* log10(ifelse(deg$p_val_adj > 0, deg$p_val_adj, correct_value))
deg$pct <- ifelse(deg$avg_log2FC > 0, deg$pct.1, deg$pct.2)
deg$gene_change <- ifelse((deg$avg_log2FC > 0.5 & deg$p_val_adj < 0.05), "Up",
                          ifelse((deg$avg_log2FC < -0.5 & deg$p_val_adj < 0.05), "Down", "Else"))
table(deg$gene_change)
epi.deg.tmp <- filter(epi.deg, p_val_adj < 0.05)
# 设置需要展示的基因名称
epi.deg.tmp <- epi.deg.tmp %>% 
  dplyr::arrange(desc(avg_log2FC))
genes <- c("CD74","HLA-DRB1","HLA-DRA","HLA-DPA1","HLA-DRB5","HLA-DPB1","HLA-DMA","HLA-DQB1","HLA-DQA1",
           "CXCL10","CD274",
           "HSPA1A","CEACAM5","PLA2G10")
deg$gene_name <- ifelse(rownames(deg) %in% genes, rownames(deg), NA)
p <- deg %>% 
  dplyr::arrange(desc(avg_log2FC)) %>%
  ggplot(aes(x = avg_log2FC, y = logp, color = gene_change)) +
  geom_point(aes(size = pct), position = jitter, alpha = 1) +  # 使用带有纵向抖动的 jitter
  geom_text_repel(aes(x = avg_log2FC, y = logp, label = gene_name),
                  color = "black", position = jitter, size = 3, fontface = "bold", 
                  max.overlaps = 50, segment.color = "#808080") +  # 设置线的颜色为灰色
  theme_bw(base_size = 15) +
  ggtitle(paste0("DEGs of DC in MSI CRC TME")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +  # log2FC = -0.3 和 0.3
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p = 0.01
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 15),
        text = element_text(face = "bold")) +
  scale_color_manual(values=c("Down"="#00468BFF", "Else"="grey", "Up"="#ED0000FF"))
p

# fig8f
genelist <- c("CD274","HLA-DRA","CXCL10","PLA2G10","HSPA1A")
# 创建一个列表来存储每个 FeaturePlot
plots <- lapply(genelist, function(gene) {
  FeaturePlot(epi.nmf, features = gene, pt.size = 1, split.by = "State") +
    theme(strip.text = element_blank())
})
plots
combined_plot <- wrap_plots(plots, ncol = 1, nrow = 5)
combined_plot

# fig8g
library(AUCell)
auc_score <- list(
  # "cGAS_STING_Score" = c("CGAS","TMEM173","TBK1","IKBKE","IRF3"),
  "IFNA_Response_Score" = dplyr::filter(H_geneset, term == "HALLMARK_INTERFERON_ALPHA_RESPONSE")[,2],
  "IFNG_Response_Score" = dplyr::filter(H_geneset, term == "HALLMARK_INTERFERON_GAMMA_RESPONSE")[,2],
  # "IFN_MHCII_Score" = sting_ISG[sting_ISG$path == "IFN_MHCII",]$gene,
  "STING_Score" = sting_score)
cells_ranking <- AUCell_buildRankings(epi.nmf@assays$RNA@data)
cells_auc <- AUCell_calcAUC(auc_score, cells_ranking, aucMaxRank = 0.1 * nrow(cells_ranking))
aa <- as.data.frame(t(cells_auc@assays@data$AUC))
tmp_info <- select(epi.nmf@meta.data, "PID","State","seurat_clusters")
tmp_info <- merge(tmp_info, aa, by = 0)
tmp_info <- pivot_longer(tmp_info, cols = c(5,6), names_to = "IFN_Response_Pathway", values_to = "IFN_Score")
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info, aes(x = State, y = IFN_Score, fill = State)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  scale_y_continuous(name = "IFN Response Score")+
  scale_x_discrete(name = "State") +
  # stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    # axis.text.x = element_blank(),  # 隐藏横轴文字
    # axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ms_colors) +
  facet_wrap(~IFN_Response_Pathway) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig8h
tmp_info <- dplyr::select(epi.nmf@meta.data, "PID", "STING_Score", "State")
compaired <- list(c("MSI","MSS"))
ggplot(data = tmp_info, aes(x = State, y = STING_Score, fill = State)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  scale_y_continuous(name = "cGAS-STING Pathway Activation Score")+
  scale_x_discrete(name = "State") +
  # stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    # axis.text.x = element_blank(),  # 隐藏横轴文字
    # axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ms_colors) +
  # facet_wrap(~STING_Score) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig8j
ggplot(data = hscore_info, aes(x = State, y = Sting, fill = State)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  scale_y_continuous(name = "STING IHC H-score")+
  scale_x_discrete(name = "State") +
  # stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    # axis.text.x = element_blank(),  # 隐藏横轴文字
    # axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ms_colors) +
  # facet_wrap(~STING_Score) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig8l
ggplot(data = hscore_info, aes(x = State, y = Pla2g10, fill = State)) +
  geom_violin(scale = "width", color = "white") +  # 小提琴图边框设置为透明
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +  # 箱线图边框透明
  scale_y_continuous(name = "Pla2g10 IHC H-score")+
  scale_x_discrete(name = "State") +
  # stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey", size = 0.5) +  # 绘制连接箱线图中点的折线
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_blank(),  # 去掉 X 轴的标题
    # axis.text.x = element_blank(),  # 隐藏横轴文字
    # axis.ticks.x = element_blank(),  # 隐藏横轴刻度
    legend.position = "none"  # 移除图例
  ) +
  scale_fill_manual(values = ms_colors) +
  # facet_wrap(~STING_Score) +
  geom_signif(
    comparisons = compaired,
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

# fig8m
med_pla <- median(ihc_res_t_pla2g10$PLA2G10, na.rm = TRUE)
ihc_res_t_pla2g10 <- ihc_res_t_pla2g10 %>%
  mutate(
    PLA2G10_group = ifelse(PLA2G10 >= med_pla, "High", "Low")
  )
ihc_res_t_pla2g10$PLA2G10_group <- factor(ihc_res_t_pla2g10$PLA2G10_group, levels = c("Low", "High"))

ggplot(ihc_res_t_pla2g10, aes(x = PLA2G10_group, y = CD8, fill = PLA2G10_group)) +
  geom_violin(aes(fill = PLA2G10_group), scale = "width", color = "white") +
  geom_boxplot(width = 0.1, position = position_identity(), outlier.alpha = 0, fill = "white", color = "black") +
  theme_classic() +
  theme(
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Low" = "#3c537a",
                               "High" = "#ee9798")) +
  geom_signif(
    comparisons = list(c("Low", "High")),
    step_increase = 0.3,
    map_signif_level = FALSE,
    test = wilcox.test
  )

lines <- readLines("/home/kongdeyang/MSI_CRC/MSIcode_fig.R")
writeLines(lines, "/home/kongdeyang/MSI_CRC/myscript.txt")
