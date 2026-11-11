library(Seurat)
args <- commandArgs(TRUE)
root <- if (length(args)) args[[1]] else "/Users/paul/WorkInProgress/dietary_intervention"
data <- file.path(root, "data", "rds")
object <- readRDS(file.path(data, "tcr_integrated.seurat.rds"))
rna <- object[["RNA"]]
layers <- Layers(rna, search = "^counts($|\\.)")
stopifnot(length(layers) > 0L)
if (length(layers) > 1L) rna <- JoinLayers(rna, layers = "counts", new = "counts")
layer <- if (length(layers) > 1L) "counts" else layers[[1]]
counts <- if (inherits(rna, "Assay")) slot(rna, layer) else LayerData(rna, layer = layer)
umap <- Embeddings(object, reduction = "umap")
rm(object, rna)
invisible(gc())
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
for (start in seq.int(1, length(counts@x), by = 1000000)) {
  values <- counts@x[seq.int(start, min(start + 999999, length(counts@x)))]
  stopifnot(all(is.finite(values)), all(values >= 0), all(values == floor(values)))
}
meta <- readRDS(file.path(data, "tcr_meta_only.rds"))
tcell <- readRDS(file.path(data, "v2a_meta_with_exh.rds"))
fields <- c("sample_id", "group", "celltype_l0", "celltype_l1", "celltype_l2")
cells <- colnames(counts)
stopifnot(!anyDuplicated(cells), all(cells %in% rownames(meta)),
          all(cells %in% rownames(umap)), all(rownames(tcell) %in% cells))
meta <- meta[cells, fields, drop = FALSE]
for (field in c("label_exh", "label_v2a")) {
  meta[[field]] <- as.character(tcell[[field]][match(cells, rownames(tcell))])
}
meta$x <- umap[cells, 1]
meta$y <- umap[cells, 2]
stopifnot(all(complete.cases(meta[, c("sample_id", "group", "x", "y")])),
          !anyDuplicated(unique(meta[, c("sample_id", "group")])$sample_id))
rm(umap, tcell)
invisible(gc())
message(nrow(meta), " cells · ", nrow(counts), " genes · ", length(unique(meta$sample_id)), " mice")
target <- Sys.getenv("EXPLORER_CACHE", unset = "explorer-data.rds")
temporary <- tempfile(tmpdir = dirname(target))
saveRDS(list(meta = meta, counts = counts), temporary, compress = FALSE)
stopifnot(file.rename(temporary, target))
