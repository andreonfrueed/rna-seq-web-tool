#!/usr/bin/env Rscript
# ============================================================
#  DESeq2 VST 后处理：VST 归一化矩阵 + 样本聚类热图 + DEG 热图
#  用法：
#    Rscript DESeq2_vst.R --counts <Raw_Counts.xlsx> --samples <samples.tsv> \
#      --gtf <gtf> --outdir <output_dir> [--top 50] [--topvar 1000] \
#      [--gene-labels 1|0]（默认 1：标注 top 显著基因）
# ============================================================

# ---------------------------------------------------------------- 统一学术配色
# Nature Publishing Group 风格：低饱和、印刷友好、色盲可辨。
PAL_UP    <- "#E64B35"   # 上调（暖红）
PAL_DOWN  <- "#4DBBD5"   # 下调（青蓝）
PAL_NS    <- "#B0B0B0"   # 不显著（浅灰）
PAL_GROUP <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
# 样本聚类热图配色（绿-黑-红）
CLUSTER_LOW  <- "#00A651"
CLUSTER_MID  <- "#000000"
CLUSTER_HIGH <- "#FF1A1A"
# DEG 热图配色（蓝-白-红）
DEG_LOW  <- "#4393C3"
DEG_MID  <- "#F7F7F7"
DEG_HIGH <- "#E64B35"

# ---------------------------------------------------------------- GTF 基因符号映射
parse_gtf_symbols <- function(gtf_path) {
  # 从 GTF 离线解析 gene_id -> gene_name，返回 data.frame(gene, symbol)。
  # 用 grep+awk+cut 只提取 gene 行（约数万行）的第 9 列属性串，
  # 避免把整个 GTF（含数百万 exon/transcript 行）读进内存。
  if (is.null(gtf_path) || !file.exists(gtf_path)) {
    return(data.frame(gene = character(0), symbol = character(0),
                      stringsAsFactors = FALSE))
  }
  cmd <- sprintf("grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"' | cut -f9",
                 shQuote(gtf_path))
  attrs <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(attrs) == 0L) {
    return(data.frame(gene = character(0), symbol = character(0),
                      stringsAsFactors = FALSE))
  }
  gid <- sub('.*gene_id "([^"]+)".*', '\\1', attrs)
  gn  <- sub('.*gene_name "([^"]+)".*', '\\1', attrs)
  keep <- gid != attrs & gn != attrs  # sub 未匹配时返回原串，据此剔除
  df <- data.frame(gene = gid[keep], symbol = gn[keep], stringsAsFactors = FALSE)
  df <- df[!duplicated(df$gene), , drop = FALSE]
  df
}

# 手写 95% 置信椭圆（二维正态协方差椭圆），不依赖 ellipse 包
cov_ellipse <- function(x, y, level = 0.95, n = 100) {
  m <- c(mean(x), mean(y))
  S <- stats::cov(cbind(x, y))
  e <- eigen(S)
  r <- sqrt(stats::qchisq(level, df = 2))
  t <- seq(0, 2 * pi, length.out = n)
  pts <- t(m + r * (e$vectors %*% diag(sqrt(e$values)) %*% rbind(cos(t), sin(t))))
  list(x = pts[, 1], y = pts[, 2])
}

parse_args <- function(argv) {
  args <- list(counts = NULL, samples = NULL, gtf = NULL, outdir = NULL,
               symbol_map = NULL, top = 50L, topvar = 1000L,
               gene_labels = TRUE)
  known <- c("--counts", "--samples", "--gtf", "--outdir", "--symbol-map",
             "--top", "--topvar", "--gene-labels")
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[[i]]
    if (!flag %in% known) {
      stop(paste0("未知参数: ", flag), call. = FALSE)
    }
    if (i + 1L > length(argv)) {
      stop(paste0("参数缺少值: ", flag), call. = FALSE)
    }
    name <- gsub("-", "_", sub("^--", "", flag))
    value <- argv[[i + 1L]]
    if (name %in% c("top", "topvar")) {
      value <- suppressWarnings(as.integer(value))
      if (is.na(value) || value <= 0L) {
        stop(paste0("参数 ", flag, " 必须是正整数"), call. = FALSE)
      }
    }
    if (name == "gene_labels") {
      value <- tolower(as.character(value))
      if (value %in% c("1", "true")) {
        value <- TRUE
      } else if (value %in% c("0", "false")) {
        value <- FALSE
      } else {
        stop(paste0("参数 ", flag, " 必须是 1/0/true/false"), call. = FALSE)
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

read_symbol_map <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(data.frame(gene = character(0), symbol = character(0),
                      stringsAsFactors = FALSE))
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("gene", "symbol") %in% names(df))) {
    stop("symbol-map 需要 gene 和 symbol 两列", call. = FALSE)
  }
  df
}

symbols_for <- function(genes, symbol_map) {
  out <- symbol_map$symbol[match(genes, symbol_map$gene)]
  out[is.na(out)] <- genes[is.na(out)]
  # 重复 symbol 加后缀，保证热图行名唯一
  dup <- duplicated(out)
  if (any(dup)) {
    out[dup] <- paste0(out[dup], "_", genes[dup])
  }
  out
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

render_deg_heatmap <- function(vst_mat, outdir, top, symbol_map) {
  diff_dir <- file.path(outdir, "4.Differential_Expression")
  deg_xlsx <- file.path(diff_dir, "Filtered_DEGs.xlsx")
  tryCatch({
    sig_genes <- character(0)
    # BUG-19：DEG 热图的显著基因来源，首选 R 自己生成的 DESeq2 CSV——
    # 差异表(diff_genes)、火山图、MA 图都读它，是 deseq2 引擎的权威来源。
    # Filtered_DEGs.xlsx 是 pyseqrna 上游产物，仅作无 CSV 时的回退。
    # 且 CSV 分支显著标准与 run_deseq_results / 火山图完全一致
    # （padj<0.05 且 |log2FC|>=1），避免同一份数据选出两套基因。
    csv_files <- list.files(diff_dir, pattern = "^DESeq2_.*_vs_.*\\.csv$",
                            full.names = TRUE)
    if (length(csv_files) > 0L) {
      for (f in csv_files) {
        deg_df <- utils::read.csv(f, stringsAsFactors = FALSE)
        if (!all(c("padj", "log2FoldChange") %in% names(deg_df))) next
        keep <- !is.na(deg_df$padj) & !is.na(deg_df$log2FoldChange) &
          deg_df$padj < 0.05 & abs(deg_df$log2FoldChange) >= 1
        sig_genes <- unique(c(sig_genes, as.character(deg_df$Gene[keep])))
      }
    } else if (file.exists(deg_xlsx)) {
      # 回退：无 DESeq2 CSV 时读 pyseqrna 的 Filtered_DEGs.xlsx；
      # p 值列 padj 优先、fdr 次之，消除对列顺序的依赖。
      sheets <- readxl::excel_sheets(deg_xlsx)
      for (sheet in sheets) {
        raw_deg <- readxl::read_excel(deg_xlsx, sheet = sheet)
        deg_df <- as.data.frame(raw_deg, stringsAsFactors = FALSE)
        names(deg_df) <- names(raw_deg)
        if (!("Gene" %in% names(deg_df))) next
        pvalue_cols <- names(deg_df)[grepl("^(fdr|padj)", tolower(names(deg_df)))]
        if (length(pvalue_cols) == 0L) next
        pvalue_cols <- c(pvalue_cols[grepl("^padj", tolower(pvalue_cols))],
                         pvalue_cols[!grepl("^padj", tolower(pvalue_cols))])
        pvalues <- deg_df[[pvalue_cols[1L]]]
        keep <- !is.na(pvalues) & pvalues < 0.05
        sig_genes <- unique(c(sig_genes, as.character(deg_df$Gene[keep])))
      }
    }
    if (length(sig_genes) == 0L) {
      warning("差异表中没有显著基因，跳过 DEG 热图", call. = FALSE)
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
    deg_pdf <- file.path(heat_dir, "DEG_heatmap_vst.pdf")
    heat_mat <- vst_mat[matched_genes, , drop = FALSE]
    rownames(heat_mat) <- symbols_for(rownames(heat_mat), symbol_map)
    draw_deg_heat <- function(filename) {
      pheatmap::pheatmap(
        heat_mat,
        scale = "row",
        cluster_cols = FALSE,  # BUG-20：样本列按 samples.tsv 顺序固定排列，不按表达相似度打乱
        main = NA,
        color = grDevices::colorRampPalette(
          c(DEG_LOW, DEG_MID, DEG_HIGH))(100),
        border_color = NA,
        fontsize = 10,
        filename = filename,
        width = 10, height = 8, res = 300
      )
    }
    draw_deg_heat(deg_png)  # 网页预览
    draw_deg_heat(deg_pdf)  # 期刊矢量投稿
  }, error = function(e) {
    warning(paste0("DEG 热图生成失败，已跳过: ", conditionMessage(e)), call. = FALSE)
  })
  invisible(NULL)
}

run_deseq_results <- function(dds, condition_levels, symbol_map, outdir) {
  diff_dir <- file.path(outdir, "4.Differential_Expression")
  gene_dir <- file.path(diff_dir, "diff_genes")
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
  # 清理 pyseqrna 内置 deseq2 写的旧 diff_genes（命名体系不同：LPS-C vs C-LPS），
  # 避免同一比较两套文件并存——R 是 deseq2 引擎的权威差异来源
  old_txt <- list.files(gene_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(old_txt) > 0L) {
    file.remove(old_txt)
  }
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  for (i in seq_len(length(condition_levels) - 1L)) {
    for (j in (i + 1L):length(condition_levels)) {
      c1 <- condition_levels[[i]]
      c2 <- condition_levels[[j]]
      res <- DESeq2::results(dds, contrast = c("condition", c1, c2))
      df <- as.data.frame(res)
      df$Gene <- rownames(df)
      df$Symbol <- symbols_for(df$Gene, symbol_map)
      df <- df[, c("Gene", "Symbol", "baseMean", "log2FoldChange",
                   "lfcSE", "stat", "pvalue", "padj")]
      rownames(df) <- NULL
      write.csv(df, file.path(diff_dir, paste0("DESeq2_", c1, "_vs_", c2, ".csv")),
                row.names = FALSE, fileEncoding = "UTF-8")
      sig <- !is.na(df$padj) & df$padj < 0.05 & !is.na(df$log2FoldChange)
      up <- df$Gene[sig & df$log2FoldChange >= 1]
      down <- df$Gene[sig & df$log2FoldChange <= -1]
      writeLines(up, file.path(gene_dir, paste0(c1, "-", c2, "_up.txt")))
      writeLines(down, file.path(gene_dir, paste0(c1, "-", c2, "_down.txt")))
      writeLines(c(up, down), file.path(gene_dir, paste0(c1, "-", c2, ".txt")))
      cat("  DESeq2", c1, "vs", c2, ": up", length(up), "down", length(down), "\n")
    }
  }
}

# ============================================================
# 出图双格式 helper：PNG(300dpi 网页预览) + PDF(矢量投稿)。
# 期刊投稿要求矢量图（PDF/SVG，缩放不糊、文字可选）；
# 同一绘图代码画两遍，分别落到两种设备。
# ============================================================
save_figure <- function(png_path, pdf_path, width, height, res = 300, plot_fn) {
  grDevices::png(png_path, width = width, height = height, units = "in", res = res)
  plot_fn()
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = width, height = height)
  plot_fn()
  grDevices::dev.off()
  invisible(NULL)
}

render_pca <- function(vst_mat, col_data, outdir) {
  pca_dir <- file.path(outdir, "5.Visualization", "Sample_Plots")
  dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
  row_vars <- apply(vst_mat, 1L, stats::var)
  n <- min(1000L, nrow(vst_mat))
  top_genes <- names(sort(row_vars, decreasing = TRUE))[seq_len(n)]
  pca <- stats::prcomp(t(vst_mat[top_genes, , drop = FALSE]),
                       center = TRUE, scale. = FALSE)
  percent <- round(100 * summary(pca)$importance[2L, ], 1)
  conds <- as.character(col_data$condition)
  lv <- unique(conds)
  colors <- PAL_GROUP[as.integer(factor(conds, levels = lv))]
  pca_png <- file.path(pca_dir, "All_Samples_PCA_vst.png")
  pca_pdf <- file.path(pca_dir, "All_Samples_PCA_vst.pdf")
  save_figure(pca_png, pca_pdf, 7.2, 6, 300, function() {
    graphics::par(mar = c(4.2, 4.2, 1.2, 4.5))
    graphics::plot(pca$x[, 1L], pca$x[, 2L], col = colors, pch = 16, cex = 1.5,
                   xlab = paste0("PC1 (", percent[[1L]], "%)"),
                   ylab = paste0("PC2 (", percent[[2L]], "%)"),
                   cex.lab = 1.3, cex.axis = 1.1, las = 1)
    # 95% 置信椭圆（按组，手写协方差椭圆，不依赖额外包）
    for (g in seq_along(lv)) {
      idx <- conds == lv[g]
      if (sum(idx) >= 2) {
        xy <- pca$x[idx, 1:2, drop = FALSE]
        tryCatch({
          ell <- cov_ellipse(xy[, 1], xy[, 2])
          graphics::lines(ell$x, ell$y, col = PAL_GROUP[g], lwd = 1.2, lty = 2)
        }, error = function(e) NULL)
      }
    }
    # 样本标签（点上方小字号，展示样本名，避免与点/椭圆重叠）
    graphics::text(pca$x[, 1L], pca$x[, 2L], labels = rownames(pca$x),
                   cex = 0.62, pos = 3, offset = 0.45, col = "#333333")
    graphics::legend("topright", legend = lv,
                     col = PAL_GROUP[seq_along(lv)],
                     pch = 16, pt.cex = 1.4, bty = "n", inset = -0.14, xpd = NA,
                     title = "Group", title.adj = 0)
  })
}

render_volcanoes <- function(outdir, gene_labels) {
  diff_dir <- file.path(outdir, "4.Differential_Expression")
  if (!dir.exists(diff_dir)) {
    return(invisible(NULL))
  }
  csv_files <- list.files(diff_dir, pattern = "^DESeq2_.*_vs_.*\\.csv$",
                          full.names = TRUE)
  for (f in csv_files) {
    tryCatch({
      deg_df <- utils::read.csv(f, stringsAsFactors = FALSE)
      if (!all(c("log2FoldChange", "padj") %in% names(deg_df))) {
        warning(paste0("火山图缺少统计列，已跳过: ", basename(f)), call. = FALSE)
        next
      }
      keep <- !is.na(deg_df$log2FoldChange) & !is.na(deg_df$padj) &
        is.finite(deg_df$log2FoldChange) & is.finite(deg_df$padj)
      plot_df <- deg_df[keep, , drop = FALSE]
      if (nrow(plot_df) == 0L) {
        warning(paste0("火山图没有有效数据，已跳过: ", basename(f)), call. = FALSE)
        next
      }
      lfc <- plot_df$log2FoldChange
      neg_log_padj <- -log10(pmax(plot_df$padj, 1e-300))
      sig <- plot_df$padj < 0.05 & abs(lfc) >= 1
      up <- sig & lfc >= 1
      down <- sig & lfc <= -1

      base <- basename(f)
      parts <- regmatches(base, regexec("^DESeq2_(.*)_vs_(.*)\\.csv$", base))[[1L]]
      if (length(parts) != 3L) {
        warning(paste0("火山图文件名无法解析，已跳过: ", base), call. = FALSE)
        next
      }
      c1 <- parts[[2L]]
      c2 <- parts[[3L]]
      volcano_dir <- file.path(outdir, "5.Visualization", "Volcano")
      dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
      volcano_png <- file.path(volcano_dir,
                               paste0(c1, "_vs_", c2, "_volcano.png"))
      volcano_pdf <- file.path(volcano_dir,
                               paste0(c1, "_vs_", c2, "_volcano.pdf"))
      label_df <- NULL
      if (isTRUE(gene_labels) && "Symbol" %in% names(plot_df) && any(sig)) {
        score <- neg_log_padj + abs(lfc)
        top_idx <- order(score, decreasing = TRUE)
        top_idx <- top_idx[sig[top_idx]]
        top_idx <- top_idx[seq_len(min(10L, length(top_idx)))]
        if (length(top_idx) > 0L) {
          yspan <- diff(range(neg_log_padj))
          gap <- max(0.07 * yspan, 0.6)
          px <- py <- tx <- ty <- sym <- side <- character(0)
          for (side_name in c("up", "down")) {
            idx <- if (side_name == "up") {
              top_idx[lfc[top_idx] > 0]
            } else {
              top_idx[lfc[top_idx] < 0]
            }
            if (length(idx) == 0L) next
            idx <- idx[order(neg_log_padj[idx], decreasing = TRUE)]
            anchor_y <- max(neg_log_padj[idx]) + 0.06 * yspan
            pos_y <- anchor_y - (seq_along(idx) - 1L) * gap
            pos_x <- if (side_name == "up") max(lfc) + 0.9 else min(lfc) - 0.9
            px <- c(px, lfc[idx])
            py <- c(py, neg_log_padj[idx])
            tx <- c(tx, rep(pos_x, length(idx)))
            ty <- c(ty, pos_y)
            sym <- c(sym, as.character(plot_df$Symbol[idx]))
            side <- c(side, rep(side_name, length(idx)))
          }
          label_df <- data.frame(px = as.numeric(px), py = as.numeric(py),
                                 tx = as.numeric(tx), ty = as.numeric(ty),
                                 sym = sym, side = side,
                                 stringsAsFactors = FALSE)
        }
      }
      save_figure(volcano_png, volcano_pdf, 8, 6.5, 300, function() {
        graphics::par(mar = c(7.2, 5.0, 1.2, 1.2))
        xlim <- range(lfc)
        ylim <- range(neg_log_padj)
        if (!is.null(label_df) && nrow(label_df) > 0L) {
          xlim <- range(c(xlim, label_df$tx))
          ylim <- range(c(ylim, label_df$ty))
          xlim <- xlim + c(-1, 1) * (0.15 * diff(xlim) + 0.5)
          ylim <- ylim + c(-1, 1) * (0.10 * diff(ylim))
        }
        # 分层绘制：先画不显著（浅灰垫底），再画上/下调（显著点置顶不被盖住）
        graphics::plot(lfc[!sig], neg_log_padj[!sig], col = PAL_NS, pch = 16,
                       cex = 0.45, xlim = xlim, ylim = ylim,
                       xlab = "log2 Fold Change",
                       ylab = "-log10(adjusted p-value)",
                       cex.lab = 1.3, cex.axis = 1.1, las = 1)
        graphics::points(lfc[up], neg_log_padj[up], col = PAL_UP, pch = 16, cex = 0.55)
        graphics::points(lfc[down], neg_log_padj[down], col = PAL_DOWN, pch = 16, cex = 0.55)
        graphics::abline(v = c(-1, 1), col = "gray60", lty = 2)
        graphics::abline(h = -log10(0.05), col = "gray60", lty = 2)
        graphics::legend("bottom", inset = c(0, -0.22), xpd = NA, horiz = TRUE,
                         legend = c(sprintf("Up (%d)", sum(up)),
                                    sprintf("Down (%d)", sum(down)),
                                    sprintf("Not significant (%d)", sum(!sig))),
                         col = c(PAL_UP, PAL_DOWN, PAL_NS),
                         pch = 16, bty = "n")
        if (!is.null(label_df) && nrow(label_df) > 0L) {
          graphics::segments(label_df$px, label_df$py,
                             label_df$tx, label_df$ty,
                             col = "gray45", lwd = 0.7)
          left_idx <- label_df$side == "down"
          draw_label <- function(idx, adj) {
            if (any(idx)) {
              # 白描边：先画白色加粗底字，再画深色细字，避免文字与点/线重叠难读
              graphics::text(label_df$tx[idx], label_df$ty[idx],
                             labels = label_df$sym[idx],
                             cex = 0.62, col = "white", font = 2, adj = adj)
              graphics::text(label_df$tx[idx], label_df$ty[idx],
                             labels = label_df$sym[idx],
                             cex = 0.62, col = "gray20", font = 1, adj = adj)
            }
          }
          draw_label(left_idx, c(1, 0.5))
          draw_label(!left_idx, c(0, 0.5))
        }
      })
    }, error = function(e) {
      if (grDevices::dev.cur() > 1L) {
        grDevices::dev.off()
      }
      warning(paste0("火山图生成失败，已跳过: ", conditionMessage(e)), call. = FALSE)
    })
  }
  invisible(NULL)
}

render_ma_plots <- function(outdir) {
  diff_dir <- file.path(outdir, "4.Differential_Expression")
  if (!dir.exists(diff_dir)) {
    return(invisible(NULL))
  }
  csv_files <- list.files(diff_dir, pattern = "^DESeq2_.*_vs_.*\\.csv$",
                          full.names = TRUE)
  for (f in csv_files) {
    tryCatch({
      deg_df <- utils::read.csv(f, stringsAsFactors = FALSE)
      if (!all(c("baseMean", "log2FoldChange", "padj") %in% names(deg_df))) {
        next
      }
      keep <- !is.na(deg_df$baseMean) & !is.na(deg_df$log2FoldChange) &
        is.finite(deg_df$baseMean) & is.finite(deg_df$log2FoldChange) &
        deg_df$baseMean > 0
      plot_df <- deg_df[keep, , drop = FALSE]
      if (nrow(plot_df) == 0L) {
        next
      }
      padj <- plot_df$padj
      padj[is.na(padj)] <- 1
      sig <- padj < 0.05 & abs(plot_df$log2FoldChange) >= 1
      base <- basename(f)
      parts <- regmatches(base, regexec("^DESeq2_(.*)_vs_(.*)\\.csv$", base))[[1L]]
      if (length(parts) != 3L) {
        next
      }
      c1 <- parts[[2L]]
      c2 <- parts[[3L]]
      ma_dir <- file.path(outdir, "5.Visualization", "MA_plots")
      dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)
      ma_png <- file.path(ma_dir, paste0(c1, "_vs_", c2, "_MA.png"))
      ma_pdf <- file.path(ma_dir, paste0(c1, "_vs_", c2, "_MA.pdf"))
      save_figure(ma_png, ma_pdf, 8, 6.5, 300, function() {
        graphics::par(mar = c(5.0, 5.0, 1.2, 1.2))
        x <- log2(plot_df$baseMean)
        y <- plot_df$log2FoldChange
        # 分层绘制：不显著垫底，显著点置顶
        up <- sig & plot_df$log2FoldChange > 0
        down <- sig & plot_df$log2FoldChange < 0
        nonsig <- !sig
        graphics::plot(x[nonsig], y[nonsig], col = PAL_NS, pch = 16, cex = 0.4,
                       xlab = "log2(mean expression)", ylab = "log2 Fold Change",
                       cex.lab = 1.3, cex.axis = 1.1, las = 1)
        graphics::points(x[up], y[up], col = PAL_UP, pch = 16, cex = 0.5)
        graphics::points(x[down], y[down], col = PAL_DOWN, pch = 16, cex = 0.5)
        graphics::abline(h = 0, col = "gray50", lty = 2)
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) > 20) {
          smooth <- stats::lowess(x[ok], y[ok], f = 0.3)
          graphics::lines(smooth, col = "#2F618C", lwd = 1.8)
        }
        graphics::legend("topright", inset = 0.02,
                         legend = c(sprintf("Up (%d)", sum(up)),
                                    sprintf("Down (%d)", sum(down)),
                                    "Not significant"),
                         col = c(PAL_UP, PAL_DOWN, PAL_NS),
                         pch = 16, bty = "n", cex = 0.9)
      })
    }, error = function(e) {
      if (grDevices::dev.cur() > 1L) {
        grDevices::dev.off()
      }
      warning(paste0("MA 图生成失败，已跳过: ", conditionMessage(e)), call. = FALSE)
    })
  }
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

  identifiers <- as.character(
    samples$samples$Identifier[match(keep_cols, samples$sample_names)])
  # 与网页端"处理÷对照"策略一致：按样本表首现顺序取逆序作为比较方向
  condition_levels <- rev(unique(identifiers))
  col_data <- data.frame(
    condition = factor(identifiers, levels = condition_levels),
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
  # BUG-24：样本列按 condition_levels 重排（首位即对照/参考组，处理组在后），
  # 组内按样本名排序；col_data 同步重排，保证 PCA 颜色与样本列一一对应。
  # 热图/VST 矩阵列序因此统一为「对照组 → 处理组」，符合阅读直觉与「对照 vs 处理」命名。
  col_order <- rownames(col_data)[order(
    match(as.character(col_data$condition), condition_levels),
    rownames(col_data)
  )]
  vst_mat <- vst_mat[, col_order, drop = FALSE]
  col_data <- col_data[col_order, , drop = FALSE]

  norm_dir <- file.path(args$outdir, "4.Normalization")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  vst_df <- as.data.frame(vst_mat, check.names = FALSE, stringsAsFactors = FALSE)
  vst_df <- cbind(Gene = rownames(vst_mat), vst_df, stringsAsFactors = FALSE)
  symbol_map <- read_symbol_map(args$symbol_map)
  # BUG-FIX：run_pipeline.sh 未传 --symbol-map 时 Symbol 列会退化成 Ensembl ID
  # （如 ENSMUSG00000000001），火山图/热图只能显示基因 ID 而非基因符号。
  # 回退为直接从 GTF 解析 gene_id -> gene_name，基因名从此可读、可投稿。
  if (nrow(symbol_map) == 0L) {
    symbol_map <- parse_gtf_symbols(args$gtf)
  }
  vst_df <- cbind(Gene = vst_df$Gene,
                  Symbol = symbols_for(vst_df$Gene, symbol_map),
                  vst_df[, -1L, drop = FALSE],
                  stringsAsFactors = FALSE)
  rownames(vst_df) <- NULL
  utils::write.csv(vst_df, file.path(norm_dir, "VST_normalized_counts.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")

  run_deseq_results(dds, condition_levels, symbol_map, args$outdir)
  render_volcanoes(args$outdir, args$gene_labels)
  render_ma_plots(args$outdir)
  render_pca(vst_mat, col_data, args$outdir)

  if (nrow(vst_mat) >= 2L) {
    row_vars <- apply(vst_mat, 1L, stats::var)
    top_var_n <- min(args$topvar, nrow(vst_mat))
    top_var_genes <- names(sort(row_vars, decreasing = TRUE))[seq_len(top_var_n)]
    cluster_dir <- file.path(args$outdir, "5.Clustering")
    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
    cluster_png <- file.path(cluster_dir, "sample_clustering_vst_heatmap.png")
    cluster_pdf <- file.path(cluster_dir, "sample_clustering_vst_heatmap.pdf")
    heat_mat <- vst_mat[top_var_genes, , drop = FALSE]
    rownames(heat_mat) <- symbols_for(rownames(heat_mat), symbol_map)
    draw_cluster_heat <- function(filename) {
      pheatmap::pheatmap(
        heat_mat,
        scale = "row",
        clustering_method = "average",
        cluster_cols = FALSE,  # BUG-23：样本列统一固定顺序，与 DEG 热图一致，不再按相似度打乱
        show_rownames = FALSE,
        main = NA,
        color = grDevices::colorRampPalette(
          c(CLUSTER_LOW, CLUSTER_MID, CLUSTER_HIGH))(100),
        border_color = NA,
        fontsize = 10,
        filename = filename,
        width = 10, height = 8, res = 300
      )
    }
    draw_cluster_heat(cluster_png)  # 网页预览
    draw_cluster_heat(cluster_pdf)  # 期刊矢量投稿
  } else {
    warning("VST 矩阵基因太少，跳过样本聚类热图", call. = FALSE)
  }

  render_deg_heatmap(vst_mat, args$outdir, args$top, symbol_map)
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
