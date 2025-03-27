#' Read files from bigscape directory
#'
#' @param big_scape_dir Directory of bigscape output such as 2024-05-14_16-00-37_hybrids_glocal
#' @param cutoff default NULL, set this when you have multiple cutoffs
#' @param reassign_GCF default TRUE. When some BGCs were assigned into different GCFs ,reassign them into the GCF with the largest number of BGCs.
#'
#' @return big_scape object
#' @export
#'
read_big_scape_dir <- function(big_scape_dir = "2024-05-14_16-00-37_hybrids_glocal",
                               cutoff = NULL, reassign_GCF = TRUE) {
  all_bgc <- readr::read_delim(file.path(big_scape_dir, "/Network_Annotations_Full.tsv"), show_col_types = FALSE)

  # 导入各个文件夹
  types_dir <- list.dirs(big_scape_dir, recursive = FALSE)
  types <- basename(types_dir)
  list.files(types_dir[1], pattern = ".network") -> tmp_net
  gsub(".*_c(.*).network", "\\1", tmp_net) -> tmp_net
  message("Found ", length(tmp_net), " networks with different cutoff:")
  message(paste0(tmp_net, collapse = ", "))
  if (!is.null(cutoff)) {
    if (cutoff[1] %in% tmp_net) {
      tmp_net <- cutoff
    } else {
      stop("No ", cutoff[1], " network found in the directory")
    }
  }
  message("Use the cutoff: ", tmp_net[1])
  cutoff <- tmp_net[1]
  if (length(tmp_net) > 1) message("Set `cutoff = ` for different read")

  # 读取网络
  dabiao("Reading networks")
  networks <- list()
  for (i in seq_along(types_dir)) {
    type_dir <- types_dir[i]
    edge_df <- readr::read_delim(paste0(types_dir[i], "/", types[i], "_c", cutoff, ".network"),
      delim = "\t", show_col_types = FALSE
    )
    networks[[types[i]]] <- edge_df[, 1:3]
  }

  # 读取GCF
  dabiao("Reading GCFs")
  GCFs <- list()
  for (i in seq_along(types_dir)) {
    type_dir <- types_dir[i]
    GCF_df <- readr::read_delim(paste0(types_dir[i], "/", types[i], "_clustering_c", cutoff, ".tsv"),
      delim = "\t", show_col_types = FALSE
    )
    colnames(GCF_df) <- c("BGC", "Family_Number")
    GCF_df <- data.frame(Type = types[i], GCF_df)
    GCFs[[types[i]]] <- GCF_df
  }

  # 整理GCF
  GCF_df <- do.call(rbind, GCFs) %>% data.frame(row.names = NULL)
  id_num <- nchar(max(GCF_df$Family_Number))
  GCF_df$GCF <- paste0("GCF", sprintf(paste0("%0", id_num, "d"), GCF_df$Family_Number))

  dplyr::distinct(GCF_df[, c("BGC", "GCF")]) %>% dplyr::count(GCF) -> GCF_count
  GCF_count <- GCF_count %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(Group = ifelse(n > 1, "GCF", "Singleton"))
  GCF_df <- dplyr::left_join(GCF_df, GCF_count)

  BGC_GCF_df <- dplyr::distinct(GCF_df[, c("BGC", "GCF", "n", "Group")]) %>%
    dplyr::arrange(-n)
  GCFs <- split(BGC_GCF_df$BGC, BGC_GCF_df$GCF)
  if (any(duplicated(BGC_GCF_df$BGC))) {
    message(
      "Some BGCs such as '", BGC_GCF_df[duplicated(BGC_GCF_df$BGC), ][1, 1],
      "' were assigned into different GCFs"
    )
  }
  if (reassign_GCF) {
    # 优先认为一个BGC属于更大的GCF
    message("Reassign them into the GCF with the largest number of BGCs")
    BGC_GCF_df <- dplyr::distinct(BGC_GCF_df, BGC, .keep_all = TRUE)
  }
  all_bgc <- dplyr::left_join(all_bgc[, 1:7], BGC_GCF_df)

  if (any(grepl("^BGC\\d{7}$", all_bgc$BGC))) {
    message("Some BGC names are in the format of BGCXXXXXXX (from MIBiG database) ")
    with_mibig <- TRUE
  } else {
    with_mibig <- FALSE
  }
  if (with_mibig) {
    all_bgc$Source <- ifelse(grepl("^BGC\\d{7}$", all_bgc$BGC), "MIBiG", "User")
  }

  big_scape_res <- list(all_bgc = all_bgc, networks = networks, GCFs = GCFs)
  class(big_scape_res) <- "big_scape"
  attributes(big_scape_res)$cutoff <- cutoff
  attributes(big_scape_res)$with_mibig <- with_mibig
  big_scape_res
}

#' Transform the network to metanet object
#'
#' @param big_scape_res big_scape object from `read_big_scape_dir()`
#' @param class default "NRPS", set this to the class you want to extract
#'
#' @return metanet object
#' @export
#'
trans_net <- function(big_scape_res, class = "NRPS") {
  stopifnot(inherits(big_scape_res, "big_scape"))
  lib_ps("MetaNet", "igraph", library = FALSE)

  net <- MetaNet::c_net_from_edgelist(big_scape_res$networks[[class]])

  bgc_anno <- big_scape_res$all_bgc %>%
    dplyr::select(BGC, `Product Prediction`, `BiG-SCAPE class`) %>%
    dplyr::distinct(BGC, .keep_all = TRUE) %>%
    tibble::column_to_rownames("BGC")
  igraph::V(net)$degree <- igraph::degree(net)
  net <- MetaNet::c_net_set(net, bgc_anno, vertex_class = "Product Prediction", vertex_size = "degree")
  igraph::E(net)$e_type <- rep(
    paste0("Distance<", attributes(big_scape_res)$cutoff),
    length(igraph::E(net))
  )
  net
}

#' Plot BGC class
#'
#' @param big_scape_res big_scape object from `read_big_scape_dir()`
#' @param mode 1~2, 1 for doughnut plot, 2 for sankey plot
#' @param rm_mibig logical, remove MIBiG BGCs from the plot
#' @param ... additional parameters for `gghuan` or `my_sankey`, such as `topN=10`
#'
#' @return ggplot
#' @export
#'
plot_BGC_class <- function(big_scape_res, mode = 1, rm_mibig = FALSE, ...) {
  stopifnot(inherits(big_scape_res, "big_scape"))
  all_bgc <- big_scape_res$all_bgc
  if (rm_mibig) {
    all_bgc <- dplyr::filter(all_bgc, Source != "MIBiG")
  }

  if (mode == 1) {
    NP_class <- dplyr::count(all_bgc, `BiG-SCAPE class`, name = "count") %>% dplyr::arrange(-count)
    p <- do.call(pcutils::gghuan, update_param(list(tab = NP_class, topN = 15, percentage = TRUE), list(...))) +
      annotate("text", 0, 0, label = paste0(format(nrow(all_bgc), big.mark = ","), " BGCs")) +
      scale_fill_pc("col2")
  } else if (mode == 2) {
    dplyr::count(all_bgc, `BiG-SCAPE class`, `Product Prediction`) %>%
      dplyr::arrange(`BiG-SCAPE class`, -n) -> NP_class2
    p <- do.call(
      pcutils::my_sankey,
      update_param(
        list(
          test = NP_class2, topN = 15,
          D3_params = list(width = 400, height = 400, numberFormat = "")
        ),
        list(...)
      )
    )
  }
  return(p)
}

#' Read BiG-SLiCE output dir
#'
#' @param big_slice_dir path to the BiG-SLiCE output directory
#'
#' @return big_slice object
#' @export
#'
read_big_slice_dir <- function(big_slice_dir) {
  filepath <- normalizePath(big_slice_dir)
  expected_structure <- c("app", "result", "requirements.txt", "start_server.sh", "result/data.db")
  if (!check_directory_structure(filepath, expected_structure = expected_structure)) {
    check_directory_structure(filepath, expected_structure = expected_structure, verbose = TRUE)
    stop("The directory structure is not correct. Please check the directory again.")
  }
  lib_ps("RSQLite", "DBI", library = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(filepath, "result/data.db"))

  # all_bgc
  bgc <- DBI::dbReadTable(con, "bgc")
  gcf <- DBI::dbReadTable(con, "gcf_membership")

  colnames(bgc)[1] <- "bgc_id"
  gcf$gcf_id <- paste0("GCF_", gcf$gcf_id)
  colnames(gcf)[3] <- "distance"
  gcfs <- dplyr::left_join(bgc[, 1:6], gcf[, 1:3])

  dplyr::count(gcfs, gcf_id) -> GCF_count
  GCF_count <- GCF_count %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(Group = ifelse(n > 1, "GCF", "Singleton"))
  all_bgc <- dplyr::left_join(gcfs, GCF_count)
  all_bgc

  GCFs <- split(all_bgc$bgc_id, all_bgc$gcf_id)

  reports <- NULL
  if (file.exists(file.path(filepath, "reports/reports.db"))) {
    con2 <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(filepath, "reports/reports.db"))
    reports <- DBI::dbReadTable(con2, "reports")
  }

  big_slice_res <- list(all_bgc = all_bgc, GCFs = GCFs, reports = reports)
  class(big_slice_res) <- "big_slice"
  attributes(big_slice_res)$tables <- DBI::dbListTables(con)
  attributes(big_slice_res)$bgc_number <- nrow(all_bgc)
  attributes(big_slice_res)$filepath <- filepath

  DBI::dbDisconnect(con)
  DBI::dbDisconnect(con2)
  big_slice_res
}

#' Open BiG-SLiCE output website
#'
#' @param big_slice_res big_slice object
#' @param port port
#'
#' @return No value
#'
open_big_slice_website <- function(big_slice_res, port = 1234) {
  stopifnot(inherits(big_slice_res, "big_slice"))
  filepath <- attributes(big_slice_res)$filepath
  system(paste("bash", file.path(filepath, "start_server.sh"), port), wait = FALSE)
  utils::browseURL(paste0("http://10.197.87.125:", port))
}


#' Get big-slice database
#'
#' @param big_slice_res big_slice object
#' @param table default "bgc"
#'
#' @return data.frame
#' @export
#'
get_big_slice_db <- function(big_slice_res, table = "bgc") {
  stopifnot(inherits(big_slice_res, "big_slice"))
  lib_ps("RSQLite", "DBI", library = FALSE)
  table <- match.arg(table, attributes(big_slice_res)$tables)
  filepath <- attributes(big_slice_res)$filepath
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(filepath, "result/data.db"))
  df <- DBI::dbReadTable(con, table)
  DBI::dbDisconnect(con)
  df
}

#' Get report dataframe from big_slice_res
#'
#' @param big_slice_res big_slice object
#' @param report report name
#' @param distance the threshold to define a bgc in a gcf, default 0.4.
#'
#' @return data.frame
#' @export
get_report_df <- function(big_slice_res, report, distance = 0.4) {
  stopifnot(inherits(big_slice_res, "big_slice"))
  lib_ps("RSQLite", "DBI", library = FALSE)
  if (is.null(big_slice_res$reports)) {
    stop("No report available.")
  }
  if (!report %in% big_slice_res$reports$name) {
    stop(paste0("Report '", report, "' not found."))
  }
  report_id <- reports[reports$name == report, "id"]

  filepath <- attributes(big_slice_res)$filepath
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(filepath, "reports", report_id, "data.db"))

  # all_bgc
  bgc <- DBI::dbReadTable(con, "bgc")
  gcf <- DBI::dbReadTable(con, "gcf_membership")

  colnames(bgc)[1] <- "bgc_id"
  gcf$gcf_id <- paste0("GCF_", gcf$gcf_id)
  colnames(gcf)[3] <- "distance"
  gcfs <- dplyr::left_join(bgc[, 1:6], gcf[, 1:3])

  gcfs$in_gcf <- ifelse(gcfs$distance < distance, TRUE, FALSE)
  gcfs
}
