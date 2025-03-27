#' Read json file from antismash directory
#'
#' @param BGC_json_file antismash output json file
#'
#' @return antismash object
#' @export
#'
read_antismash_json <- function(BGC_json_file) {
  pcutils::lib_ps("jsonlite", library = FALSE)
  if (dir.exists(BGC_json_file)) {
    BGC_json_file <- list.files(BGC_json_file, ".json")
  }
  BGC_json <- jsonlite::read_json(BGC_json_file)

  BGC_json$genome_name <- gsub(".json", "", basename(BGC_json_file))
  find_BGC <- vapply(BGC_json$records, \(i)i$id %in% names(BGC_json$timings), logical(1))
  BGC_json$records <- BGC_json$records[find_BGC]
  names(BGC_json$records) <- vapply(BGC_json$records, \(i)i$id, character(1))
  class(BGC_json) <- "antismash"
  BGC_json
}


parse_location <- function(x = "[0:8952](+)") {
  start <- as.numeric(gsub(x = x, pattern = "\\[(\\d+):(\\d+)\\].*", replacement = "\\1"))
  end <- as.numeric(gsub(x = x, pattern = "\\[(\\d+):(\\d+)\\].*", replacement = "\\2"))
  if (grepl("\\(", x)) {
    direct <- gsub(x = x, pattern = ".*\\((.)\\)", replacement = "\\1")
  } else {
    direct <- "+"
  }
  return(c(start = start, end = end, direct = ifelse(direct == "+", 1, -1)))
}

get_region_df <- function(single_contig) {
  region_site <- vapply(single_contig[["features"]], \(i)i[["type"]] == "region", logical(1)) %>% which()
  single_contig_region_res <- list()
  for (pos in region_site) {
    single_region <- single_contig[["features"]][[pos]]
    single_contig_region_res[[as.character(pos)]] <- data.frame(
      id = single_contig$id,
      contig = single_contig$description,
      region_number = single_region$qualifiers$region_number[[1]],
      on_contig_edge = single_region$qualifiers$contig_edge[[1]],
      length = diff(parse_location(single_region$location)[1:2]) %>% unname(),
      product = single_contig[["features"]][[pos]]$qualifiers$product %>% paste0(collapse = "; ")
    )
  }
  do.call(rbind, single_contig_region_res)
}

#' Get information from antismash output json file
#'
#' @param BGC_json antismash object from `read_antismash_json`
#'
#' @return data.frame
#' @export
#'
get_BGCs_from_BGC_json <- function(BGC_json) {
  stopifnot(inherits(BGC_json, "antismash"))

  BGC_res <- BGC_json$records
  genome_name <- BGC_json$genome_name
  if (length(BGC_res) > 0) {
    BGC_product_df <- list()
    for (BGC_i in seq_along(BGC_res)) {
      BGC_product_df[[BGC_i]] <- get_region_df(BGC_res[[BGC_i]])
    }
    BGC_product_df <- do.call(rbind, BGC_product_df)
    return(data.frame(genome_name = genome_name, BGC_product_df, row.names = NULL))
  } else {
    warning("No BGCs in ", genome_name)
    return(data.frame())
  }
}

get_features_df <- function(single_contig) {
  features_df <- data.frame(
    id = rep(single_contig$id, length(single_contig$features)),
    description = rep(single_contig$description, length(single_contig$features))
  )
  lapply(single_contig$features, \(i)data.frame(type = i$type, location = i$location)) %>%
    do.call(rbind, .) -> tmp_df
  features_df <- cbind(features_df, tmp_df, lapply(tmp_df$location, parse_location) %>% do.call(rbind, .))
  features_df$start <- features_df$start + 1
  features_df
}

get_cds_df <- function(single_contig) {
  get_features_df(single_contig) -> tmp_df
  cds_df <- dplyr::filter(tmp_df, type == "CDS")
  single_contig_cds <- single_contig$features[tmp_df$type == "CDS"]
  cds_df$locus_tag <- names(single_contig_cds) <- vapply(single_contig_cds, \(i)i$qualifiers$locus_tag[[1]], character(1))
  lapply(single_contig_cds, \(i)i$qualifiers$gene_kind[[1]]) -> tmpls
  tmpls[unlist(lapply(tmpls, is.null))] <- "other"
  cds_df$gene_kind <- unlist(tmpls)
  lapply(single_contig_cds, \(i)paste0(unlist(i$qualifiers$gene_functions), collapse = "; ")) -> tmpls
  tmpls[unlist(lapply(tmpls, is.null))] <- NA
  cds_df$gene_functions <- unlist(tmpls)
  cds_df
}

BGC_gene_kind_col <- setNames(
  c("#8dd3c7", "#F8CC00", "#bebada", "#fb8072", "#80b1d3", "grey"),
  c(
    "biosynthetic-additional", "transport", "regulatory",
    "biosynthetic", "resistance", "other"
  )
)


#' Plot a BGC
#'
#' @param BGC_json antismash object from `read_antismash_json`
#' @param region_id region_id, default is NULL
#' @param show_locus_tag logical, default is TRUE
#' @param show_gene_functions logical, default is FALSE
#'
#' @return ggplot
#' @export
#'
plot_BGC <- function(BGC_json, region_id = NULL, show_locus_tag = TRUE,
                     show_gene_functions = FALSE) {
  stopifnot(inherits(BGC_json, "antismash"))

  if (is.null(region_id)) {
    message("region_id is NULL, use the first region as default")
    region_id <- names(BGC_json$records)[1]
  }
  single_contig <- BGC_json$records[[region_id]]
  get_region_df(single_contig) -> region_df

  get_cds_df(single_contig) -> cds_df
  lib_ps("gggenes", library = FALSE)

  p <- ggplot2::ggplot(data = cds_df, aes(xmin = start, xmax = end, y = id))
  p <- p + gggenes::geom_gene_arrow(
    mapping = aes(fill = gene_kind, forward = (direct != -1))
  ) + scale_fill_manual(values = BGC_gene_kind_col)

  if (show_locus_tag) {
    p <- p + gggenes::geom_gene_label(
      mapping = aes(label = locus_tag), min.size = 0
    )
  }
  if (show_gene_functions) {
    lib_ps("ggrepel", "stringr", library = FALSE)
    p <- p +
      ggrepel::geom_label_repel(
        aes(
          x = c(.data$start + .data$end) / 2, fill = gene_kind,
          label = stringr::str_wrap(gene_functions, 30)
        ),
        box.padding = 0.6, size = 3, nudge_y = c(1, -1), segment.curvature = 0.01,
        label.r = 0, show.legend = FALSE
      )
  }

  sub_title <- single_contig$description
  p <- p + coord_fixed(diff(range(c(cds_df$start, cds_df$start))) / 8)
  p <- p + gggenes::theme_genes() +
    labs(
      y = NULL, title = paste0(
        "Genome: ", BGC_json$genome_name,
        "; Region: ", region_id,
        "; Product: ", paste0(region_df$product, collapse = ", ")
      ),
      subtitle = sub_title
    ) +
    theme(
      legend.position = "top", axis.line.y = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
  p
}


get_knowncluster_from_BGC_json <- function(BGC_json_file) {
  pcutils::lib_ps("jsonlite", library = FALSE)
  BGC_json <- jsonlite::read_json(BGC_json_file)

  genome_name <- gsub(".json", "", basename(BGC_json_file))
  # 找到region的contigs有哪些：
  # names(BGC_json$timings)

  find_BGC <- vapply(BGC_json$records, \(i)i$id %in% names(BGC_json$timings), logical(1))
  BGC_res <- BGC_json$records[find_BGC]

  if (length(BGC_res) > 0) {
    knowncluster_df <- list()
    for (BGC_i in seq_along(BGC_res)) {
      knowncluster_df[[BGC_i]] <- from_contig_get_knowncluster(BGC_res[[BGC_i]])
    }
    knowncluster_df <- do.call(rbind, knowncluster_df)

    if (!is.null(knowncluster_df)) {
      if (!"similarity" %in% colnames(knowncluster_df)) {
        knowncluster_df$similarity <- floor(100 * knowncluster_df$n_hit_protein / knowncluster_df$n_protein)
      }
      return(data.frame(genome_name = genome_name, knowncluster_df))
    } else {
      warning("No Knowncluster in ", genome_name)
      return(data.frame())
    }
  } else {
    warning("No BGCs in ", genome_name)
    return(data.frame())
  }
}

from_contig_get_knowncluster <- function(single_contig) {
  knowncluster <- single_contig$modules$antismash.modules.clusterblast$knowncluster$results
  single_contig_region_res <- list()
  for (i in seq_along(knowncluster)) {
    tmp_res <- knowncluster[[i]]
    if (tmp_res$total_hits > 0) {
      single_contig_region_res[[i]] <- data.frame(
        id = single_contig$id,
        contig = single_contig$description,
        region_number = tmp_res$region_number,
        lapply(tmp_res$ranking, \(i){
          c(i[[1]][-c(3, 6)],
            n_protein = length(i[[1]]$proteins),
            n_hit_protein = lapply(i[[2]]$pairings, \(i)i[[3]]$name) %>% unique() %>% length(),
            i[[2]][-6]
          ) %>%
            as.data.frame()
        }) %>% do.call(rbind, .)
      )
    }
  }
  knowncluster_df <- do.call(rbind, single_contig_region_res)
}
