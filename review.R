read_reviews <- function(path) {
  if (file.exists(path)) {
    return(read.csv(path, colClasses = "character", na.strings = NULL, check.names = FALSE))
  }
  data.frame(cluster = character(), label = character(), evidence = character(),
             status = character(), changed_at = character())
}

save_review <- function(path, cluster, label, evidence, status) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(paste0(path, ".lock"), timeout = 3000)
  if (is.null(lock)) stop("다른 창에서 주석을 저장 중입니다. 다시 시도하세요.")
  on.exit(filelock::unlock(lock), add = TRUE)
  history <- rbind(read_reviews(path), data.frame(
    cluster = cluster, label = label, evidence = evidence, status = status,
    changed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3 UTC", tz = "UTC")
  ))
  temp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(temp), add = TRUE)
  write.csv(history, temp, row.names = FALSE, fileEncoding = "UTF-8")
  if (!file.rename(temp, path)) stop("주석 파일을 저장할 수 없습니다.")
  history
}

reviewed_labels <- function(meta, history) {
  current <- history[!duplicated(history$cluster, fromLast = TRUE), ]
  current <- current[current$status == "applied", ]
  replacement <- current$label[match(meta$cluster_id, current$cluster)]
  labels <- as.character(meta$celltype_l1)
  labels[!is.na(replacement)] <- replacement[!is.na(replacement)]
  labels
}

cluster_markers <- function(counts, in_cluster, gene_totals) {
  block <- counts[, in_cluster, drop = FALSE]
  inside <- Matrix::rowSums(block)
  outside <- gene_totals - inside
  if (sum(inside) == 0 || sum(outside) == 0) stop("비교할 RNA 카운트가 충분하지 않습니다.")
  detection <- tabulate(block@i[block@x > 0] + 1L, nbins = nrow(block)) / ncol(block)
  log_ratio <- log2((1 + 1e6 * inside / sum(inside)) / (1 + 1e6 * outside / sum(outside)))
  keep <- detection >= 0.1 & log_ratio > 0
  markers <- data.frame(gene = rownames(counts)[keep], log2_ratio = log_ratio[keep],
                        detected_pct = 100 * detection[keep])
  markers[order(-markers$log2_ratio, -markers$detected_pct, markers$gene), ]
}
