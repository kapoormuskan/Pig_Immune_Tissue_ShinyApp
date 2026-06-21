#.libPaths("/work/ABG/MUSKAN/ShinyTissues/micromamba/envs/seurat5.3.0/lib/R/library") Sys.setenv(RSTUDIO_WHICH_R = "/work/ABG/MUSKAN/ShinyTissues/micromamba/envs/seurat5.3.0/bin/R")

#loading lib paths
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)
library(glmGamPoi)

# Auto-detect output directory, directly saves to app's data/directory
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
output_dir <- normalizePath(file.path(script_dir, "../data"), mustWork = FALSE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


#load reference(example: BM) .rds file and query(example: SP) .rds datasets
ref <- readRDS(file.path(script_dir, "../data/PigImmuneMultiTissue_BoneMarrow_CellTypeAnnotated_HBBhiRemoved.rds"))
query <- readRDS(file.path(script_dir, "../data/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.rds"))

# Use SCTransform or perform SCTransform on both datasets
#SCT normalizes library soze effects and is used for anchor based integration, this wilol retian all genes in SCT assay 
ref <- SCTransform(ref, return.only.var.genes = FALSE, verbose = TRUE)
query <- SCTransform(query, return.only.var.genes = FALSE, verbose = TRUE)
#using default assay as SCT 
DefaultAssay(ref) <- "SCT"
DefaultAssay(query) <- "SCT"
# Set cell types in reference (in this case bone marrow)
ref$celltypes <- Idents(ref)

# Split query(spleen here) by sample ID
#label transfer is performed per sample to avoid batch effects
query.list <- SplitObject(query, split.by = "ID")
# Normalize and run PCA for each sample
query.list <- lapply(query.list, function(obj) {
  obj <- SCTransform(obj, return.only.var.genes = FALSE, verbose = TRUE)
  obj <- RunPCA(obj, npcs = 50)
  return(obj)
})
CellTypePredictions <- list()
MappingScores <- list()

#Anchor finding and label tramsfer. Performing mapping for each query spleen sample

#FindTransferAnchors: Identifies mutual nearest neigbors between reference and query in CCA space, then these anchors are used to transfer labels and compute mappign scores
for (i in seq_along(query.list)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query.list[[i]],
    reduction = "cca",
    dims = 1:30,
    normalization.method = "SCT",
    recompute.residuals = FALSE
  )
#TransferData: uses the anchors to predict cell labels for each query cell and retyrns a score for every refrecen cell type plus a predicted ID (highest score) 
  predictions <- TransferData(
    anchorset = anchors,
    refdata = list(cell_type = ref$celltypes),
    dims = 1:30,
    weight.reduction = "cca"
  )
#MappingScore: measure of how well each query cell is placed in the rederebce space >.5 is good (range is 0-1)
  map_scores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors = query.list[[i]]@neighbors[["query.list_ref.nn"]],
    query.weights = Tool(query.list[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(query.list[[i]]),
    ref.embeddings = Embeddings(ref),
    nn.method = "annoy"
  )
  CellTypePredictions[[i]] <- predictions
  MappingScores[[i]] <- data.frame(MappingScores = map_scores)
}

# we combine prediction scores + mapping scores 
predictions.df <- do.call(rbind, CellTypePredictions) %>% as.data.frame()
predictions.df$CellBarcodes <- rownames(predictions.df)
mapping.df <- do.call(rbind, MappingScores)
mapping.df$CellBarcodes <- rownames(mapping.df)
colnames(mapping.df)[1] <- "MappingScores"
#Merging predictions & mapping scores by cell barcodes (musta have same cells in same order)
merged_meta <- merge(predictions.df, mapping.df, by = "CellBarcodes")
rownames(merged_meta) <- merged_meta$CellBarcodes

#run PCA again on SCT assay
query <- RunPCA(query, assay = "SCT", npcs = 50)
query <- RunUMAP(query, dims = 1:30, reduction = "pca")

#now add metadata to seurat with default assay RNA
DefaultAssay(query) <- "RNA"
query <- AddMetaData(query, metadata = merged_meta)


#saving results
write_xlsx(predictions.df, file.path(output_dir, "BoneMarrowToSpleen_CellTypePredictions.xlsx"))
write_xlsx(mapping.df, file.path(output_dir, "BoneMarrowToSpleen_MappingScores.xlsx"))
saveRDS(query, file = file.path(output_dir, "BoneMarrowToSpleen_Annotated.rds"))

#visualiztion of mapping scores - QC
FeaturePlot(query, features = "MappingScores", reduction = "umap", min.cutoff = 0.4) + scale_color_gradientn(colors = c('white', 'yellow', 'orange', 'red'))
#UMAP for predicted cell types
DimPlot(query, group.by = "predicted.id", reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Predicted Cell Types (from Reference)")
#violin plot of mapping scores plotted per predicted cell type
VlnPlot(query, features = "MappingScores", group.by = "predicted.id") + ggtitle("Mapping Score by Predicted Cell Type")
#how confidently it matches to predicted cell type?
VlnPlot(query, features = "prediction.score.max", group.by = "predicted.id") + ggtitle("Prediction Score (max) by Cell Type")

#plotting
library(patchwork) #predcition scores 
score_cols <- grep("^prediction\\.score\\.", colnames(query@meta.data), value = TRUE)
score_cols_filtered <- score_cols[sapply(score_cols, function(x) max(query[[x]][, 1]) > 0.1)]
#plotting featyre plot
plot_list <- lapply(score_cols_filtered, function(feature) {FeaturePlot(query, features = feature, reduction = "umap") + scale_color_gradientn(colors = c("grey90", "red")) + NoAxes() + ggtitle(gsub("prediction.score.", "", feature)) + theme(plot.title = element_text(size = 8, hjust = 0.5))  # smaller title
})
wrap_plots(plot_list, ncol = 5)


####################################################
#libararies loading
library(dplyr)
library(ggplot2)
library(scales)

#creating dataframes for dot plot: Size= % of original cells predicted as the cel;l type and color is the average mapping score
df <- data.frame(
CellType = query$CellType,
PredictedID = query$predicted.id,
MappingScore = query$MappingScores
)
#plotting counts and percentage
count_df <- df %>%
group_by(CellType, PredictedID) %>%
summarise(Count = n(), .groups = "drop")
pct_df <- count_df %>%
group_by(CellType) %>%
mutate(Percent = 100 * Count / sum(Count)) %>%
ungroup()

#average mapping scores
score_df <- df %>%
group_by(CellType, PredictedID) %>%
summarise(MappingScore = mean(MappingScore), .groups = "drop")

#% +mapping score
dot_df <- merge(pct_df, score_df, by = c("CellType", "PredictedID"))

ggplot(dot_df, aes(x = PredictedID, y = CellType)) +
geom_point(aes(size = Percent, fill = MappingScore), shape = 21, color = "black") +
scale_size_continuous(name = "Percent (%)", range = c(1, 7)) +
scale_fill_gradientn(colors = c("yellow", "orange", "red", "darkred"),
limits = c(0.5, 1), oob = squish,
name = "Avg Mapping Score") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
axis.text.y = element_text(size = 10),
axis.title = element_text(size = 12)) +
labs(title = "Predicted vs Annotated Cell Types",
x = "Predicted Cell Type",
y = "Original Annotation")
