#' @title Get MiBIG database files
#'
#' @param download_dir where to save the download .tar.gz file?
#' @param file mibig_json_.*tar.gz from https://dl.secondarymetabolites.org/mibig/
#'
#' @return No value
#' @export
get_mibig_db <- function(download_dir, file = NULL) {
  pcutils::lib_ps("jsonlite", library = FALSE)
  if (is.null(file)) {
    if (missing(download_dir)) {
      stop("Please set download_dir!\ndownload_dir must be specified if file is not specified")
    }
    if (!dir.exists(download_dir)) {
      dir.create(download_dir)
      if (!dir.exists(download_dir)) stop("Can not create download_dir!")
    }
    url <- get_mibig_versions()
    file <- file.path(download_dir, basename(url))
    pcutils::download2(url, file)
  }

  pack_dir <- tools::R_user_dir("BGCkit")
  if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

  if (!grepl("mibig_json_.*tar.gz", file)) {
    stop("file name must be like 'mibig_json_3.1.tar.gz'")
  }

  mibig_version <- gsub(".*(\\d\\.\\d)\\.tar\\.gz", "\\1", file)
  pcutils::dabiao("Loading ", mibig_version, " mibig_json")

  utils::untar(tarfile = file, exdir = dirname(file))

  dir <- gsub(".tar.gz", "", file)
  json_list <- list.files(path = dir, pattern = "BGC.*\\.json")
  pcutils::dabiao("Found ", length(json_list), " BGCs\n")
  lapply(json_list, \(js)jsonlite::read_json(paste0(dir, "/", js))) -> mibig_db
  names(mibig_db) <- gsub(".json", "", json_list)
  attributes(mibig_db)$version <- mibig_version
  attributes(mibig_db)$download_time <- Sys.time()

  class(mibig_db) <- "mibig_db"

  save(mibig_db, file = paste0(pack_dir, "/mibig_db.rda"))
  pcutils::dabiao(paste0("Update done at ", Sys.time()))
}

#' Load the mibig_db
#'
#' @param verbose logical
#'
#' @export
#' @return mibig_db
load_mibig_db <- function(verbose = TRUE) {
  prefix <- "mibig_db"
  new_file <- file.path(tools::R_user_dir("BGCkit"), paste0(prefix, ".rda"))
  envir <- environment()

  if (file.exists(new_file)) {
    load(new_file, envir = envir)
  } else {
    message("Not find ", prefix, ", please run `get_", prefix, "()` first!")
    return(invisible())
  }
  res <- get(prefix, envir = envir)
  if (verbose) {
    pcutils::dabiao("load ", prefix)
    if (!is.null(attributes(res)$"download_time")) {
      pcutils::dabiao(paste0(prefix, " download time: ", attributes(res)$"download_time"))
      message("If you want to update ", prefix, ", use `get_", prefix, "()`")
    }
  }
  return(res)
}

#' @title Get MiBIG database version
#'
#' @return Url
get_mibig_versions <- function() {
  lib_ps("rvest", library = FALSE)
  url <- "https://dl.secondarymetabolites.org/mibig/"
  page <- rvest::read_html(url)
  page %>%
    rvest::html_elements("a") %>%
    rvest::html_attr("href") -> links
  links[grepl("mibig_json.*\\.gz$", links)] -> versions
  message("Finded versions:\n", paste0(versions, collapse = "\n"))
  return(paste0(url, versions[length(versions)]))
}

#' @title Convert MiBIG database to data frame
#'
#' @param mibig_db mibig_db object from `load_mibig_db()`
#'
#' @return mibig_db_df
#' @export
mibig_db2df <- function(mibig_db) {
  stopifnot(inherits(mibig_db, "mibig_db"))
  mibig_version <- attributes(mibig_db)$version
  message("MiBIG databse version: ", mibig_version)

  mibig_cluster_df <- list()
  if (numeric_version(mibig_version) == numeric_version("3.1")) {
    for (js in names(mibig_db)) {
      mibig_json <- mibig_db[[js]]
      mibig_cluster_df[[js]] <-
        data.frame(mibig_json$cluster[c("mibig_accession", "minimal", "ncbi_tax_id", "organism_name", "status")],
          publications = paste0(mibig_json$cluster$publications, collapse = "; ")
        )
      # 不一定有chem_struct，mol_mass等
      tmp_comp <-
        lapply(mibig_json$cluster$compounds, \(i){
          data.frame(
            compound = i$compound,
            chem_struct = ifelse(is.null(i$chem_struct), NA, i$chem_struct),
            mol_mass = ifelse(is.null(i$mol_mass), NA, i$mol_mass),
            molecular_formula = ifelse(is.null(i$molecular_formula), NA, i$molecular_formula),
            database_id = paste0(unlist(i$database_id), collapse = "; "),
            chem_acts = paste0(unlist(i$chem_acts), collapse = "; "),
            chem_moieties = paste0(unlist(i$chem_moieties), collapse = "; "),
            chem_targets = paste0(unlist(i$chem_targets), collapse = "; ")
          )
        }) %>% do.call(rbind, .)
      mibig_cluster_df[[js]] <- cbind(mibig_cluster_df[[js]], tmp_comp)
    }
    mibig_cluster_df <- data.frame(do.call(rbind, mibig_cluster_df), row.names = NULL)
  } else if (numeric_version(mibig_version) > numeric_version("3.1")) {
    for (js in names(mibig_db)) {
      mibig_json <- mibig_db[[js]]
      mibig_cluster_df[[js]] <-
        data.frame(mibig_json[c("accession", "quality", "completeness", "status")],
          ncbi_tax_id = mibig_json$taxonomy$ncbiTaxId,
          organism_name = mibig_json$taxonomy$name,
          publications = paste0(mibig_json$legacy_references, collapse = "; ")
        )
      # compounds可能为空
      if (length(mibig_json$compounds) == 0) {
        tmp_comp <- data.frame(
          compound = NA, chem_struct = NA, mol_mass = NA, molecular_formula = NA,
          database_id = NA, chem_acts = NA, chem_moieties = NA
        )
      } else {
        tmp_comp <-
          lapply(mibig_json$compounds, \(i){
            data.frame(
              compound = i$name,
              chem_struct = ifelse(is.null(i$structure), NA, i$structure),
              mol_mass = ifelse(is.null(i$mass), NA, i$mass),
              molecular_formula = ifelse(is.null(i$formula), NA, i$formula),
              database_id = paste0(unlist(i$databaseIds), collapse = "; "),
              chem_acts = paste0(unlist(lapply(i$bioactivities, \(i)i[["name"]])), collapse = "; "),
              chem_moieties = paste0(unlist(i$moieties), collapse = "; ")
            )
          }) %>% do.call(rbind, .)
      }
      mibig_cluster_df[[js]] <- cbind(mibig_cluster_df[[js]], tmp_comp)
    }
    mibig_cluster_df <- data.frame(do.call(rbind, mibig_cluster_df), row.names = NULL)
  }

  attributes(mibig_cluster_df)$version <- mibig_version
  mibig_cluster_df
}
