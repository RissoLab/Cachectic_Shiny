# seu <- readRDS("~/Downloads/multiome_cachectic/multiome_sgl_comb.rds")
# library(Seurat)
# sce <- as.SingleCellExperiment(seu)
# library(scater)
# 
# reducedDim(sce, "PCA", withDimnames=TRUE) <- seu[['pca']]@cell.embeddings
# reducedDim(sce, "UMAP", withDimnames=TRUE) <- seu[['umap']]@cell.embeddings
# # reducedDim(sce, "lsi", withDimnames=TRUE) <- seu[['lsi']]@cell.embeddings
# # reducedDim(sce, "integratedLSI", withDimnames=TRUE) <- seu[['integratedLSI']]@cell.embeddings
# # reducedDim(sce, "UMAP.atac", withDimnames=TRUE) <- seu[['UMAP.atac']]@cell.embeddings
# reducedDim(sce, "UMAP.wnn", withDimnames=TRUE) <- seu[['UMAP.wnn']]@cell.embeddings

# saveRDS(sce, file="~/Downloads/multiome_cachectic/multiome_sgl_comb_sce.rds")
library(iSEE)
sce <- readRDS("~/Downloads/multiome_cachectic/multiome_sgl_comb_sce.rds")
umap1 <- ReducedDimensionPlot( Type="UMAP.wnn", 
                    XAxis=2L, YAxis=1L, PanelWidth=6L, 
                    ColorByColumnData="Core.Type", 
                    ColorBy="Column data", VisualBoxOpen=FALSE)
umap1[["ColorByColumnData"]] <- "subclusters"

umap2 <- ReducedDimensionPlot( Type="UMAP.wnn", 
                    XAxis=2L, YAxis=1L, PanelWidth=6L, 
                    ColorByFeatureName="Core.Type",
                    ColorBy="Feature name", VisualBoxOpen=FALSE)
# umap2[["ColorByFeatureName"]] <- "Trim63
dt <- RowDataTable(PanelWidth=3L, Selected="Trim63")
umap2[["ColorByFeatureSource"]] <- "RowDataTable1"

fap <- FeatureAssayPlot(DataBoxOpen=FALSE, VisualChoices=c("Facet", "Color"), 
                PanelWidth=9L, FacetColumnBy="Column data",
                FacetRowBy="Column data",
                ColorByColumnData="Core.Type", 
                ColorBy="Column data",
                FacetColumnByColData="Core.Type", 
                FacetRowByColData="Core.Type",
                YAxisFeatureName="Core.Type"
                )

# iSEE(sce, 
#      initial = list(RowDataTable(Selected = "Trim63"),
#                     ReducedDimensionPlot(ColorBy = "Feature name", 
#                                          ColorByFeatureSource = "RowDataTable1"),
#                     FeatureAssayPlot(YAxisFeatureSource = "RowDataTable1")))

# VisualBoxOpen=TRUE, VisualChoices=c("Facet"),
#         FacetRowBy="Column data", FacetRowByColData="subclusters", 
#         FacetColumnBy="Column data", FacetColumnByColData="Core.Type")

fap[["FacetColumnByColData"]] <- "subclusters"
fap[["ColorByColumnData"]] <- "subclusters"
fap[["FacetRowByColData"]] <- "condition"
# fap[["YAxisFeatureName"]] <- "Trim63"
fap[["YAxisFeatureSource"]] <- "RowDataTable1"


# fapct <- fap
iSEE(sce, initial=list(dt, fap, umap1, umap2))
# iSEE(sce, initial=list(umap1, umap2, fap))
    # app
