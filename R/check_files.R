get_mibig <- function(download_dir, file = NULL) {
  if (is.null(file)) {
    if (missing(download_dir)) {
      stop("download_dir must be specified if file is not specified")
    }
    if (!dir.exists(download_dir)) {
      dir.create(download_dir)
    }
    url <- get_mibig_versions()
    file <- file.path(download_dir, basename(url))
    pcutils::download2(url, file)
  }

  if (!grepl("mibig_json_.*tar.gz", file)) {
    stop("file name must be like'mibig_json_3.1.tar.gz'")
  }

  utils::untar(tarfile = file, exdir = dirname(file))
  mibig_version <- gsub(".*(\\d\\.\\d)\\.tar\\.gz", "\\1", file)

  dir <- gsub(".tar.gz", "", file)
  json_list <- list.files(path = dir, pattern = "BGC.*\\.json")
  message("Finded ", length(json_list), " BGCs\n")
  lapply(json_list, \(js)jsonlite::read_json(paste0(dir, "/", js))) -> mibig_db
  names(mibig_db) <- gsub(".json", "", json_list)
  attributes(mibig_db)$version <- mibig_version
  save(mibig_db, file = paste0(dirname(file), "/mibig_db.rda"))
  message("Saved file: ", paste0(dirname(file), "/mibig_db.rda"))
}

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

mibig_db2df <- function(mibig_db) {
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
