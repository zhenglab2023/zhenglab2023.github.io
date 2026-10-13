# =============================================================================
# Seurat PBMC3k 分析
# 官方教程：https://satijalab.org/seurat/articles/pbmc3k_tutorial
# 官方源码：https://github.com/satijalab/seurat/blob/main/vignettes/pbmc3k_tutorial.Rmd
# =============================================================================

# =============================================================================
# 01. 加载 R 包，读取数据并创建 Seurat 对象
# =============================================================================

# 加载 dplyr、Seurat、patchwork 包
library(dplyr)
library(Seurat)
library(patchwork)

# Read10X 读取 data.dir 指定文件夹中的 10x 计数文件。
# data.dir 填写包含 matrix.mtx、genes.tsv、barcodes.tsv 的文件夹路径。
# 矩阵的行对应基因、列对应细胞，其中的数值为 UMI 计数。
pbmc.data <- Read10X(data.dir = "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/")

# CreateSeuratObject 根据计数矩阵创建 Seurat 对象，并保存为 pbmc。
# counts = pbmc.data：使用上一行读取的原始计数矩阵。
# project = "pbmc3k"：将对象的项目名称设为 pbmc3k。
# min.cells = 3：保留至少在 3 个细胞中计数大于 0 的基因。
# min.features = 200：保留至少检测到 200 个基因的细胞。
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

# =============================================================================
# 02. 查看计数矩阵并比较内存占用
# =============================================================================

# [行, 列] 从矩阵中取出指定位置的数据。
# c(...) 将三个基因名组合为行索引；1:30 表示第 1 到第 30 列细胞。
# 此行显示这三个基因在前 30 个细胞中的原始计数。
pbmc.data[c("CD3D","TCL1A","MS4A1"), 1:30]

# as.matrix 将 pbmc.data 转换为普通矩阵。
# object.size 计算转换后对象的内存大小。
dense.size <- object.size(as.matrix(pbmc.data))

# 显示普通矩阵占用的内存大小。
dense.size

# object.size 计算原始稀疏矩阵 pbmc.data 的内存大小。
sparse.size <- object.size(pbmc.data)

# 显示稀疏矩阵占用的内存大小。
sparse.size

dense.size / sparse.size

# =============================================================================
# 03. 计算质量控制指标并筛选细胞
# =============================================================================

# PercentageFeatureSet 计算每个细胞中指定基因的计数占总计数的百分比。
# pattern = "^MT-"：选择名称以 MT- 开头的基因；^ 表示字符串的开头。
# pbmc[["percent.mt"]]：在 pbmc 的细胞元数据中新增或更新 percent.mt 列。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# @meta.data 取出 pbmc 内保存细胞信息的表格。
# head(..., 5) 显示该表格的前 5 行，每行对应一个细胞。
head(pbmc@meta.data, 5)

# VlnPlot 为 features 指定的三个指标分别绘制小提琴图。
# nFeature_RNA 是检测到的基因数；nCount_RNA 是总 UMI 数；percent.mt 是上一步的百分比。
# c(...) 将三个指标名组合起来；ncol = 3 将三张图排成三列。
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter 根据两个指标绘制散点图，每个点对应一个细胞。
# feature1 指定横轴为总 UMI 数，feature2 指定纵轴为 percent.mt。
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

# 将总 UMI 数作为横轴、检测到的基因数作为纵轴，生成的图保存到 plot2。
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 使用 patchwork 的 + 将 plot1 和 plot2 并排组合并显示。
plot1 + plot2

# subset 按 subset 参数中的条件保留细胞，结果重新保存到 pbmc。
# nFeature_RNA > 200：检测到的基因数大于 200。
# nFeature_RNA < 2500：检测到的基因数小于 2500。
# percent.mt < 5：线粒体基因的计数比例小于 5%。
# & 表示“并且”，三个条件必须同时满足；等于这些边界值的细胞不保留。
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# =============================================================================
# 04. 归一化
# =============================================================================

# NormalizeData 对 pbmc 的表达数据进行归一化，结果重新保存到 pbmc。
# normalization.method = "LogNormalize"：每个细胞的基因计数先除以该细胞总计数，
# 再乘以 scale.factor，最后计算 log(1 + 处理后的数值)。
# scale.factor = 1e4：1e4 是科学计数法，等于 10000。
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)

# 不填写上述两个参数时，NormalizeData 使用它们的默认值。
# pbmc <- NormalizeData(pbmc)

# =============================================================================
# 05. 选择并显示高变基因
# =============================================================================

# FindVariableFeatures 在 pbmc 中选择高变基因，结果重新保存到 pbmc。
# selection.method = 'vst'：指定使用 vst 方法。
# nfeatures = 2000：选择 2000 个高变基因。
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

# VariableFeatures(pbmc) 取出已选定的高变基因名称。
# head(..., 10) 取其中前 10 个名称，<- 将它们保存到 top10。
top10 <- head(VariableFeatures(pbmc), 10)

# VariableFeaturePlot 根据 pbmc 中的高变基因结果生成图，保存为 plot1。
plot1 <- VariableFeaturePlot(pbmc)

# LabelPoints 在 plot 指定的图上添加 points 指定的基因名称标签。
# plot = plot1：使用上一行生成的图；points = top10：标记这 10 个基因。
# repel = TRUE：使标签尽量互相避开；结果保存为 plot2。
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 将未标注的 plot1 与已标注的 plot2 组合并显示。
plot1 + plot2

# =============================================================================
# 06. 标准化表达数据
# =============================================================================

# rownames(pbmc) 取出 pbmc 中全部基因的名称，保存为 all.genes。
all.genes <- rownames(pbmc)

# ScaleData 对 features 指定的基因进行中心化和标准化。
pbmc <- ScaleData(pbmc, features = all.genes)

# 移除不需要的变异源
# vars.to.regress = 'percent.mt'：在 ScaleData 中将 percent.mt 作为回归变量。
# pbmc <- ScaleData(pbmc, vars.to.regress = 'percent.mt')

# =============================================================================
# 07. PCA 与 PCA 结果展示
# =============================================================================

# RunPCA 对 pbmc 进行主成分分析，结果重新保存到 pbmc。
# VariableFeatures(object = pbmc) 取出 pbmc 的高变基因名称。
# features 指定将这些高变基因作为本次 PCA 的输入。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# pbmc[['pca']] 取出对象中名为 pca 的降维结果。
# print 显示结果；dims = 1:5 指定第 1 至第 5 个主成分。
# nfeatures = 5 指定每个主成分的正向和负向各显示 5 个基因。
print(pbmc[['pca']], dims = 1:5, nfeatures = 5)

# VizDimLoadings 展示基因在指定主成分上的载荷。
# dims = 1:2 指定前两个主成分；reduction = 'pca' 指定使用 PCA 结果。
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')

# DimPlot 根据 reduction = 'pca' 指定的 PCA 结果绘制细胞散点图。
# + NoLegend() 从这张图中移除图例。
DimPlot(pbmc, reduction = 'pca') + NoLegend()

# DimHeatmap 绘制指定主成分相关基因的表达热图。
# dims = 1 指定第 1 个主成分；cells = 500 指定展示的细胞数。
# balanced = TRUE 表示从正载荷和负载荷方向选取数量平衡的基因。
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# 与上一行使用相同设置，但 dims = 1:15 改为分别展示前 15 个主成分。
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# =============================================================================
# 08. 主成分概览与聚类
# =============================================================================

# ElbowPlot 根据 pbmc 的 PCA 结果绘图，横轴为主成分编号，默认纵轴为标准差。
ElbowPlot(pbmc)

# FindNeighbors 使用 dims 指定的主成分建立细胞邻居图。
# dims = 1:10：使用前 10 个主成分；结果重新保存到 pbmc。
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# FindClusters 根据 pbmc 中的邻居图进行聚类。
# resolution = 0.5：将聚类分辨率参数设为 0.5；结果重新保存到 pbmc。
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Idents(pbmc) 取出每个细胞当前的分群身份。
# head(..., 5) 显示前 5 个细胞的群编号。
head(Idents(pbmc), 5)

# =============================================================================
# 09. UMAP 与对象保存示例
# =============================================================================

# RunUMAP 计算 UMAP 降维结果，dims = 1:10 指定输入为前 10 个主成分。
pbmc <- RunUMAP(pbmc, dims = 1:10)

# DimPlot 绘制细胞散点图，reduction = 'umap' 指定使用刚计算的 UMAP 坐标。
DimPlot(pbmc, reduction = 'umap')

# saveRDS 将当前 pbmc 保存为 RDS 文件；file 指定保存路径和文件名。
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# =============================================================================
# 10. 查找标记基因
# =============================================================================

# FindMarkers 比较 ident.1 指定的群与其余细胞，结果保存为 cluster2.markers。
# ident.1 = 2 指定群 2；未指定 ident.2 时，比较组为其余所有群。
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)

# head 显示 cluster2.markers 结果表的前 5 行，n = 5 指定行数。
head(cluster2.markers, n = 5)

# 比较群 5 与群 0、群 3 合并后的细胞，结果保存为 cluster5.markers。
# ident.1 = 5 指定目标群；ident.2 = c(0, 3) 将群 0 和群 3 指定为比较组。
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))

# 显示 cluster5.markers 结果表的前 5 行。
head(cluster5.markers, n = 5)

# FindAllMarkers 依次比较每个群与其余所有细胞，合并结果保存为 pbmc.markers。
# only.pos = TRUE 表示只保留在目标群中表达较高的标记基因。
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# 从 pbmc.markers 表开始；%>% 将这一行的结果传给下一行函数。
pbmc.markers %>%
    # 按 cluster 列中的群编号分组，再通过 %>% 传给下一行。
    group_by(cluster) %>%
    # dplyr::filter 明确调用 dplyr 包中的 filter，只保留 avg_log2FC 大于 1 的行。
    # avg_log2FC 是目标群相对于比较组的平均表达差异字段；此处直接显示筛选结果。
    dplyr::filter(avg_log2FC > 1)

# FindMarkers 比较群 0 与其余细胞，结果保存为 cluster0.markers。
# logfc.threshold = 0.25：设置参与分析的基因的平均表达差异阈值。
# test.use = "roc"：指定使用 ROC 方法；only.pos = TRUE：只保留正向标记。
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# =============================================================================
# 11. 标记基因表达图与热图
# =============================================================================

# VlnPlot 为 features 指定的 MS4A1、CD79A 两个基因绘制小提琴图。
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# 为 NKG7、PF4 两个基因绘制小提琴图。
# slot = 'counts' 指定使用原始计数；log = TRUE 指定使用对数刻度。
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)

# FeaturePlot 在细胞降维散点图上用颜色表示指定基因的表达量。
# features = c(...) 将列出的 9 个基因指定为需要分别展示的基因。
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# 从全部群的标记结果表开始，%>% 将表传给下一行。
pbmc.markers %>%
    # 按 cluster 列分组，使后续操作分别作用于每个群。
    group_by(cluster) %>%
    # 在各群中保留 avg_log2FC 大于 1 的行。
    dplyr::filter(avg_log2FC > 1) %>%
    # 按当前行顺序取各群前 10 行；某群不足 10 行时取该群的全部行。
    slice_head(n = 10) %>%
    # ungroup 移除分组；-> top10 将整个处理结果赋值给右侧的 top10。
    # 此处 top10 替换了第 05 节保存的 10 个高变基因名称。
    ungroup() -> top10

# top10$gene 取出 top10 表中的 gene 列，作为热图需要展示的基因名称。
# DoHeatmap 生成这些基因的表达热图。
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# =============================================================================
# 12. 重命名细胞群
# =============================================================================

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

# levels(pbmc) 按顺序取得当前群的身份标签。
# names(new.cluster.ids) <- ... 将这些标签设置为名称向量各元素的名字，从而形成“原群标签 = 新细胞类型名称”的对应关系。
names(new.cluster.ids) <- levels(pbmc)

# RenameIdents 根据 new.cluster.ids 中的对应关系替换细胞身份标签。
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 使用 UMAP 坐标绘制细胞图；label = TRUE 显示各群的新名称。
# pt.size = 0.5 设置细胞点大小；+ NoLegend() 移除图例。
DimPlot(pbmc, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()

# =============================================================================
# 13. 设置图形格式、保存示例与环境信息
# =============================================================================

library(ggplot2)

# DimPlot 生成 UMAP 图；label = TRUE 显示群名称，label.size = 4.5 设置标签字号。
# xlab("UMAP 1") 和 ylab("UMAP 2") 分别设置横轴和纵轴的标题。
# <- 将后续各行共同设置完成的图保存到 plot；行尾的 + 表示继续添加图形设置。
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  # theme 设置外观；axis.title 指定坐标轴标题样式，legend.text 指定图例文字样式。
  # element_text(size = 18) 将这两类文字的字号设为 18。
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
  # guides 设置图例；colour 指定颜色图例，guide_legend 设置其显示形式。
  # override.aes = list(size = 10) 只将图例中的点大小设为 10。
  guides(colour = guide_legend(override.aes = list(size = 10)))


# filename 指定 JPG 文件路径；height = 7、width = 12 指定图像的高、宽（默认单位为英寸）。
# plot = plot 指定保存上面生成的图；quality = 50 设置 JPEG 图像质量。
# ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

# saveRDS 保存当前 pbmc；file 指定 RDS 文件路径和名称。
# saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

# sessionInfo 显示当前 R 版本、运行平台，以及已加载的 R 包及其版本。
sessionInfo()
