mouse_counts <- function(meta, selected, test, reference, minimum) {
  selected <- unique(selected)
  samples <- unique(meta[, c("sample_id", "group")])
  samples <- samples[samples$group %in% c(test, reference), , drop = FALSE]
  samples <- samples[order(samples$group, samples$sample_id), , drop = FALSE]
  samples$cells <- tabulate(match(meta$sample_id[selected], samples$sample_id),
                            nbins = nrow(samples))
  samples$included <- samples$cells >= minimum
  rownames(samples) <- NULL
  samples
}

run_de <- function(counts, meta, selected, samples, test, reference) {
  selected <- unique(selected)
  samples <- samples[samples$included, , drop = FALSE]
  group <- factor(samples$group, levels = c(reference, test))
  if (any(table(group) < 2L)) stop("포함된 마우스가 각 식이군에 최소 2마리 필요합니다.")
  membership <- match(meta$sample_id[selected], samples$sample_id)
  keep <- !is.na(membership)
  selected <- selected[keep]
  membership <- membership[keep]
  map <- Matrix::sparseMatrix(i = seq_along(selected), j = membership, x = 1,
                             dims = c(length(selected), nrow(samples)),
                             dimnames = list(NULL, samples$sample_id))
  bulk <- as.matrix(counts[, selected, drop = FALSE] %*% map)
  if (any(colSums(bulk) == 0)) stop("RNA 카운트가 0인 마우스가 포함되어 있습니다.")
  design <- model.matrix(~ group)
  y <- edgeR::DGEList(bulk, group = group)
  expressed <- edgeR::filterByExpr(y, design)
  if (sum(expressed) < 2L) stop("발현 필터를 통과한 유전자가 충분하지 않습니다.")
  y <- edgeR::calcNormFactors(y[expressed, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  table <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2L), n = Inf)$table
  data.frame(Gene = rownames(table), table, row.names = NULL)
}
