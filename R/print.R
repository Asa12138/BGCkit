#' Print mibig_db
#'
#' @param x mibig_db
#' @param ... Additional arguments
#'
#' @return No value
#' @exportS3Method
#' @method print mibig_db
#'
print.mibig_db <- function(x, ...) {
  # 获取 BGC 的数量
  bgc_count <- length(x)

  # 获取版本信息
  version_info <- attributes(x)$version

  # 输出结果
  cat("MiBIG Database Summary:\n")
  cat("Number of BGCs:", bgc_count, "\n")
  cat("Version:", version_info, "\n")
  cat("* this is a list, use `mibig_db2df()` to get a dataframe.")
}

#' Print antismash
#'
#' @param x antismash
#' @param ... Additional arguments
#'
#' @return No value
#' @exportS3Method
#' @method print antismash
#'
print.antismash <- function(x, ...) {
  # 获取 BGC 的数量
  bgc_count <- length(x$records)

  # 获取版本信息
  version <- x$version

  # 输出结果
  cat("Antismash Output Summary:\n")
  cat("Genome name:", x$genome_name, "\n")
  cat("Number of BGCs:", bgc_count, "\n")
  cat("Version:", version, "\n")
}


#' Print big_scape
#'
#' @param x big_scape
#' @param ... Additional arguments
#'
#' @return No value
#' @exportS3Method
#' @method print big_scape
#'
print.big_scape <- function(x, ...) {
  # 获取 BGC 的数量
  bgc_count <- nrow(x$all_bgc)

  # 获取版本信息
  cutoff <- attributes(x)$cutoff

  # 输出结果
  cat("BiG-scape Output Summary:\n")
  cat("Number of BGCs:", bgc_count, "\n")
  cat("Cutoff:", cutoff, "\n")
  cat("BiG-scape Classes:", paste0(names(x$network), collapse = ", "), "\n")
  cat("With MIBiG: ", attributes(x)$with_mibig)
}

#' Print big_slice
#'
#' @param x big_slice
#' @param ... Additional arguments
#'
#' @return No value
#' @exportS3Method
#' @method print big_slice
#'
print.big_slice <- function(x, ...) {
  # 输出结果
  cat("BiG-slice Output Summary:\n")
  cat("Dir:", attributes(x)$filepath, "\n")
  cat("BGCs number:", attributes(x)$bgc_number, "\n")
  if (!is.null(x$reports)) {
    cat("With", nrow(x$reports), "reports:\n", paste0(x$reports$name, collapse = "\n "))
  }
}

#' How to use some softwares to analysis bgcs
#'
#' @param step "antismash", "big-scape", "big-slice"
#'
#' @return No value
#' @export
#'
how_to_do_bgc <- function(step = c("antismash", "big-scape", "big-slice")) {
  step <- match.arg(step)

  if (step == "antismash") {
    res_text <- paste0(
      "# ca antismash_5.2.0 \n",
      "antismash genomes.fa --genefinding-tool prodigal \n",
      "  --cb-general --cb-subclusters --cb-knownclusters --output-dir MAG_BGC/ \n",
      "  -c 4 --minlength 10000\n"
    )
  }

  if (step == "big-scape") {
    res_text <- paste0(
      "~/miniconda3/envs/antismash_5.2.0/bin/python ~/biosoft/BiG-SCAPE-1.1.5/bigscape.py \n",
      "  -i gbk_files -o BiG_output2 \n",
      "  --min_bgc_size 10000 -c 16 -v --mibig\n"
    )
  }

  if (step == "big-slice") {
    res_text <- paste0(
      "# ca antismash_5.2.0 \n",
      "bigslice -i input_folder_template/ test_out/ \n"
    )
  }
  clipr::write_clip(res_text)
  message(res_text)
}
