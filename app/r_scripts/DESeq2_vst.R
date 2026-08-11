#!/usr/bin/env Rscript
# ============================================================
#  DESeq2 VST 后处理：VST 归一化矩阵 + 样本聚类热图 + DEG 热图
#  用法：
#    Rscript DESeq2_vst.R --counts <Raw_Counts.xlsx> --samples <samples.tsv> \
#      --gtf <gtf> --outdir <output_dir> [--top 50] [--topvar 1000]
# ============================================================

parse_args <- function(argv) {
  args <- list(counts = NULL, samples = NULL, gtf = NULL, outdir = NULL,
               top = 50L, topvar = 1000L)
  known <- c("--counts", "--samples", "--gtf", "--outdir", "--top", "--topvar")
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[[i]]
    if (!flag %in% known) {
      stop(paste0("未知参数: ", flag), call. = FALSE)
    }
    if (i + 1L > length(argv)) {
      stop(paste0("参数缺少值: ", flag), call. = FALSE)
    }
    name <- sub("^--", "", flag)
    value <- argv[[i + 1L]]
    if (name %in% c("top", "topvar")) {
      value <- suppressWarnings(as.integer(value))
      if (is.na(value) || value <= 0L) {
        stop(paste0("参数 ", flag, " 必须是正整数"), call. = FALSE)
      }
    }
    args[[name]] <- value
    i <- i + 2L
  }
  required <- c("counts", "samples", "gtf", "outdir")
  missing <- required[vapply(required, function(x) is.null(args[[x]]), logical(1L))]
  if (length(missing) > 0L) {
    stop(paste0("缺少必需参数: --", paste(missing, collapse = ", --")), call. = FALSE)
  }
  args
}

check_packages <- function() {
  required_pkgs <- c("readxl", "pheatmap", "DESeq2")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop(paste0(
      "缺少 R 包: ", paste(missing_pkgs, collapse = ", "),
      "。请安装: conda install -n pyseqrna -c conda-forge -c bioconda ",
      "r-readxl r-pheatmap bioconductor-deseq2"
    ), call. = FALSE)
  }
}

load_counts <- function(counts_path) {
  raw_counts <- readxl::read_excel(counts_path, sheet = 1L)
  counts_df <- as.data.frame(raw_counts, stringsAsFactors = FALSE)
  names(counts_df) <- names(raw_counts)
  if (ncol(counts_df) < 2L || names(counts_df)[1L] != "Gene") {
    stop("counts 表第一列必须是 Gene，其余列为样本", call. = FALSE)
  }
  gene_ids <- as.character(counts_df[[1L]])
  if (anyDuplicated(gene_ids) > 0L) {
    stop("counts 表的 Gene 列存在重复", call. = FALSE)
  }
  list(genes = gene_ids, sample_cols = names(counts_df)[-1L], data = counts_df)
}

load_samples <- function(samples_path) {
  samples_df <- utils::read.delim(
    samples_path, sep = "\t", stringsAsFactors = FALSE,
    check.names = FALSE, comment.char = ""
  )
  required_cols <- c("SampleName", "Replication", "Identifier", "File1", "File2")
  missing_cols <- setdiff(required_cols, names(samples_df))
  if (length(missing_cols) > 0L) {
    stop(paste0("samples.tsv 缺少列: ", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  sample_names <- as.character(samples_df$SampleName)
  if (anyDuplicated(sample_names) > 0L) {
    stop("samples.tsv 的 SampleName 存在重复", call. = FALSE)
  }
  list(samples = samples_df, sample_names = sample_names)
}

render_deg_heatmap <- function(vst_mat, outdir, top) {
  deg_xlsx <- file.path(outdir, "4.Differential_Expression", "Filtered_DEGs.xlsx")
  if (!file.exists(deg_xlsx)) {
    warning(paste0("未找到差异表，跳过 DEG 热图: ", deg_xlsx), call. = FALSE)
    return(invisible(NULL))
  }
  tryCatch({
    sheets <- readxl::excel_sheets(deg_xlsx)
    sig_genes <- character(0)
    for (sheet in sheets) {
      raw_deg <- readxl::read_excel(deg_xlsx, sheet = sheet)
      deg_df <- as.data.frame(raw_deg, stringsAsFactors = FALSE)
      names(deg_df) <- names(raw_deg)
      if (!("Gene" %in% names(deg_df))) {
        warning(paste0("差异表 sheet「", sheet, "」缺少 Gene 列，已跳过"), call. = FALSE)
        next
      }
      # 列名可能是 FDR 或 FDR(比较名) 等带后缀形式
      pvalue_cols <- names(deg_df)[grepl("^(fdr|padj)", tolower(names(deg_df)))]
      if (length(pvalue_cols) == 0L) {
        warning(paste0("差异表 sheet「", sheet, "」缺少 FDR/padj 列，已跳过"), call. = FALSE)
        next
      }
      pvalues <- deg_df[[pvalue_cols[1L]]]
      keep <- !is.na(pvalues) & pvalues < 0.05
      sig_genes <- unique(c(sig_genes, as.character(deg_df$Gene[keep])))
    }
    if (length(sig_genes) == 0L) {
      warning("差异表中没有 FDR/padj < 0.05 的显著基因，跳过 DEG 热图", call. = FALSE)
      return(invisible(NULL))
    }
    sig_genes <- utils::head(sig_genes, top)
    matched_genes <- intersect(sig_genes, rownames(vst_mat))
    if (length(matched_genes) == 0L) {
      warning("显著基因与 VST 矩阵没有交集，跳过 DEG 热图", call. = FALSE)
      return(invisible(NULL))
    }
    keep_var <- apply(vst_mat[matched_genes, , drop = FALSE], 1L, stats::var) > 0
    matched_genes <- matched_genes[keep_var]
    if (length(matched_genes) == 0L) {
      warning("显著基因在选中样本中方差全为 0，跳过 DEG 热图", call. = FALSE)
      return(invisible(NULL))
    }
    if (length(matched_genes) < length(sig_genes)) {
      warning("部分显著基因不在 VST 矩阵中，只绘制可匹配的基因", call. = FALSE)
    }
    heat_dir <- file.path(outdir, "5.Visualization", "Heatmaps")
    dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
    deg_png <- file.path(heat_dir, "DEG_heatmap_vst.png")
    grDevices::png(deg_png, width = 10, height = 8, units = "in", res = 300)
    pheatmap::pheatmap(
      vst_mat[matched_genes, , drop = FALSE],
      scale = "row"
    )
    grDevices::dev.off()
  }, error = function(e) {
    warning(paste0("DEG 热图生成失败，已跳过: ", conditionMessage(e)), call. = FALSE)
  })
  invisible(NULL)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  check_packages()

  for (path_arg in c("counts", "samples", "gtf")) {
    if (!file.exists(args[[path_arg]])) {
      stop(paste0("文件不存在: ", args[[path_arg]]), call. = FALSE)
    }
  }

  counts <- load_counts(args$counts)
  samples <- load_samples(args$samples)

  missing_samples <- setdiff(samples$sample_names, counts$sample_cols)
  if (length(missing_samples) > 0L) {
    stop(paste0("counts 表缺少样本列: ", paste(missing_samples, collapse = ", ")), call. = FALSE)
  }
  keep_cols <- intersect(counts$sample_cols, samples$sample_names)
  extra_cols <- setdiff(counts$sample_cols, samples$sample_names)
  if (length(extra_cols) > 0L) {
    warning(paste0("counts 表存在不在 samples.tsv 中的样本列，已忽略: ",
                   paste(extra_cols, collapse = ", ")), call. = FALSE)
  }

  counts_mat <- as.matrix(counts$data[, keep_cols, drop = FALSE])
  rownames(counts_mat) <- counts$genes
  storage.mode(counts_mat) <- "numeric"
  counts_mat[is.na(counts_mat)] <- 0
  counts_mat <- round(counts_mat)
  keep_rows <- rowSums(counts_mat > 0) > 0
  if (!any(keep_rows)) {
    stop("counts 表没有非零基因", call. = FALSE)
  }
  counts_mat <- counts_mat[keep_rows, , drop = FALSE]

  col_data <- data.frame(
    condition = factor(samples$samples$Identifier[match(keep_cols, samples$sample_names)]),
    row.names = keep_cols,
    stringsAsFactors = FALSE
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = col_data,
    design = ~condition
  )
  vsd <- DESeq2::vst(dds, blind = TRUE)
  vst_mat <- SummarizedExperiment::assay(vsd)

  norm_dir <- file.path(args$outdir, "4.Normalization")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  vst_df <- as.data.frame(vst_mat, check.names = FALSE, stringsAsFactors = FALSE)
  vst_df <- cbind(Gene = rownames(vst_mat), vst_df, stringsAsFactors = FALSE)
  rownames(vst_df) <- NULL
  utils::write.csv(vst_df, file.path(norm_dir, "VST_normalized_counts.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")

  if (nrow(vst_mat) >= 2L) {
    row_vars <- apply(vst_mat, 1L, stats::var)
    top_var_n <- min(args$topvar, nrow(vst_mat))
    top_var_genes <- names(sort(row_vars, decreasing = TRUE))[seq_len(top_var_n)]
    cluster_dir <- file.path(args$outdir, "5.Clustering")
    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
    cluster_png <- file.path(cluster_dir, "sample_clustering_vst_heatmap.png")
    grDevices::png(cluster_png, width = 10, height = 8, units = "in", res = 300)
    pheatmap::pheatmap(
      vst_mat[top_var_genes, , drop = FALSE],
      scale = "row",
      clustering_method = "average",
      show_rownames = FALSE
    )
    grDevices::dev.off()
  } else {
    warning("VST 矩阵基因太少，跳过样本聚类热图", call. = FALSE)
  }

  render_deg_heatmap(vst_mat, args$outdir, args$top)
  invisible(NULL)
}

status <- tryCatch({
  main()
  0L
}, error = function(e) {
  message(paste0("DESeq2 VST post-processing failed: ", conditionMessage(e)))
  1L
})
if (status == 0L) {
  cat("VST post-processing done\n")
}
quit(status = status)
