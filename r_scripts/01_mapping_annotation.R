args_all <- commandArgs(trailingOnly=FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=", args_all)]), winslash="/", mustWork=TRUE)
script_root <- dirname(script_path)
source(file.path(script_root, "lib", "common.R"))
check_dependencies()
source(file.path(script_root, "lib", "annotation_cache.R"))
source(file.path(script_root, "lib", "uniprot_api.R"))
source(file.path(script_root, "lib", "ncbi_api.R"))

parameters <- parse_arguments(commandArgs(trailingOnly=TRUE))
started <- format(Sys.time(), tz="UTC", usetz=TRUE)
output_root <- normalizePath(parameters$output, winslash="/", mustWork=FALSE)
tables_dir <- file.path(output_root, "tables")
raw_dir <- file.path(output_root, "raw", parameters$`run-id`)
uniprot_raw <- file.path(raw_dir, "uniprot")
ncbi_raw <- file.path(raw_dir, "ncbi")
snapshot_dir <- file.path(parameters$project, "scripts", "runs", paste0(parameters$`run-id`, "_mapping_annotation"))
if (dir.exists(raw_dir) || dir.exists(snapshot_dir)) stop("run_id já existe; snapshots não podem ser sobrescritos: ", parameters$`run-id`)
for (path in c(tables_dir, uniprot_raw, ncbi_raw, snapshot_dir)) dir.create(path, recursive=TRUE, showWarnings=FALSE)

file.copy(script_path, file.path(snapshot_dir, basename(script_path)), overwrite=FALSE)
dir.create(file.path(snapshot_dir, "lib"), showWarnings=FALSE)
for (path in list.files(file.path(script_root, "lib"), pattern="\\.R$", full.names=TRUE))
  file.copy(path, file.path(snapshot_dir, "lib", basename(path)), overwrite=FALSE)
write_json(parameters, file.path(snapshot_dir, "parameters.json"))
capture.output(sessionInfo(), file=file.path(snapshot_dir, "session_info.txt"))

input <- readr::read_csv(parameters$input, show_col_types=FALSE, name_repair="minimal")
expanded <- expand_identifiers(input, parameters$`id-column`, parameters$`id-type`)
unique_ids <- unique(expanded$input_id)
connection <- cache_connect(parameters$cache)
on.exit(DBI::dbDisconnect(connection), add=TRUE)
refresh <- identical(tolower(parameters$refresh), "true")

payloads <- list()
missing <- character()
for (id in unique_ids) {
  cached <- if (!refresh) cache_get(connection, "uniprot", id, parameters$`id-type`, parameters$`tax-id`) else NULL
  if (is.null(cached)) missing <- c(missing, id) else payloads[[id]] <- cached
}

uniprot_endpoint <- NA_character_
uniprot_release <- NA_character_
if (length(missing)) {
  fields_payload <- validate_uniprot_database(parameters$`id-type`)
  write_json(fields_payload, file.path(uniprot_raw, "idmapping_fields.json"))
  batches <- split(missing, ceiling(seq_along(missing) / 500))
  for (batch_index in seq_along(batches)) {
    batch <- batches[[batch_index]]
    job_id <- submit_uniprot_job(batch, parameters$`id-type`, parameters$`tax-id`)
    wait_uniprot_job(job_id)
    fetched <- fetch_uniprot_results(job_id)
    uniprot_endpoint <- fetched$endpoint
    uniprot_release <- fetched$release %||% NA_character_
    write_json(fetched$data, file.path(uniprot_raw, sprintf("batch_%03d.json", batch_index)))
    for (id in batch) {
      selected <- Filter(function(item) identical(as.character(item$from), id), fetched$data$results %||% list())
      payload <- list(results=selected, failedIds=if (id %in% unlist(fetched$data$failedIds %||% list())) id else list())
      cache_put(connection, "uniprot", id, parameters$`id-type`, parameters$`tax-id`, payload)
      payloads[[id]] <- payload
    }
  }
}
write_json(payloads, file.path(uniprot_raw, "effective_payloads.json"))

empty_candidate <- function(id) data.frame(input_id=id, uniprot_accession=NA_character_, uniprot_entry_name=NA_character_,
  uniprot_reviewed=NA, gene_symbol=NA_character_, uniprot_gene_names=NA_character_, organism_name=NA_character_,
  tax_id=NA_character_, protein_length=NA_integer_, protein_name=NA_character_, uniprot_function=NA_character_,
  uniprot_subcellular_location=NA_character_, ncbi_gene_id=NA_character_, ensembl_gene_ids=NA_character_,
  refseq_ids=NA_character_, stringsAsFactors=FALSE)

catalog_rows <- list()
for (row_index in seq_len(nrow(expanded))) {
  source <- expanded[row_index, , drop=FALSE]
  candidates <- parse_uniprot_response(payloads[[source$input_id]])
  resolution <- resolve_candidate_flags(candidates, parameters$`tax-id`)
  count <- resolution$candidate_count
  if (!count) candidates <- empty_candidate(source$input_id)
  preferred <- if (count) resolution$preferred_candidate else FALSE
  candidates <- cbind(source[rep(1, nrow(candidates)), ], candidates[, setdiff(names(candidates), "input_id"), drop=FALSE])
  candidates$candidate_count <- count
  candidates$mapping_status <- resolution$mapping_status
  candidates$preferred_candidate <- preferred
  candidates$preferred_reason <- ifelse(preferred, "only reviewed UniProtKB candidate", NA_character_)
  candidates$organism_match <- if (count) resolution$organism_match else FALSE
  catalog_rows[[length(catalog_rows)+1]] <- candidates
}
catalog <- if (length(catalog_rows)) do.call(rbind, catalog_rows) else cbind(expanded, empty_candidate(NA)[FALSE, -1])

ncbi_status <- "completed"
ncbi_rows <- data.frame()
ncbi_effective <- list()
gene_ids <- unique(unlist(strsplit(catalog$ncbi_gene_id[!is.na(catalog$ncbi_gene_id) & nzchar(catalog$ncbi_gene_id)], ";", fixed=TRUE)))
fetch_ncbi_cached <- function(values, mode, prefix) {
  rows <- list(); missing_values <- character()
  for (value in values) {
    cached <- if (!refresh) cache_get(connection, "ncbi", value, mode, parameters$`tax-id`) else NULL
    if (is.null(cached)) missing_values <- c(missing_values, value) else {
      ncbi_effective[[paste(mode, value, sep=":")]] <<- cached
      parsed <- parse_ncbi_response(cached)
      if (nrow(parsed)) rows[[length(rows)+1]] <- parsed
    }
  }
  if (length(missing_values)) {
    batches <- split(missing_values, ceiling(seq_along(missing_values) / 100))
    for (i in seq_along(batches)) {
      fetched <- fetch_ncbi_batch(batches[[i]], mode, parameters$`tax-id`)
      write_json(fetched$data, file.path(ncbi_raw, sprintf("%s_batch_%03d.json", prefix, i)))
      reports <- fetched$data$reports %||% list()
      for (value in batches[[i]]) {
        selected <- Filter(function(item) {
          gene <- item$gene %||% item
          if (mode == "id") identical(as.character(gene$geneId %||% ""), value)
          else identical(toupper(as.character(gene$symbol %||% "")), toupper(value))
        }, reports)
        payload <- list(reports=selected)
        ncbi_effective[[paste(mode, value, sep=":")]] <<- payload
        cache_put(connection, "ncbi", value, mode, parameters$`tax-id`, payload)
        parsed <- parse_ncbi_response(payload)
        if (nrow(parsed)) rows[[length(rows)+1]] <- parsed
      }
    }
  }
  if (length(rows)) unique(do.call(rbind, rows)) else data.frame()
}
tryCatch({
  if (length(gene_ids)) {
    ncbi_rows <- fetch_ncbi_cached(gene_ids, "id", "gene_id")
  }
  remaining_symbols <- unique(catalog$gene_symbol[!is.na(catalog$gene_symbol) & nzchar(catalog$gene_symbol)])
  if (length(remaining_symbols)) {
    symbol_rows <- fetch_ncbi_cached(remaining_symbols, "symbol", "symbol")
    if (nrow(symbol_rows)) ncbi_rows <- unique(rbind(ncbi_rows, symbol_rows))
  }
}, error=function(error) {
  ncbi_status <<- paste("annotation unavailable:", conditionMessage(error))
  writeLines(conditionMessage(error), file.path(ncbi_raw, "error.txt"))
})
write_json(ncbi_effective, file.path(ncbi_raw, "effective_payloads.json"))

ncbi_columns <- c("ncbi_official_symbol","ncbi_gene_description","ncbi_gene_type","synonyms","hgnc_identifier",
  "ncbi_ensembl_gene_ids","swissprot_accessions","ncbi_summary","ncbi_summary_source","ncbi_summary_date",
  "ncbi_tax_id","ncbi_taxonomy_name")
for (column in ncbi_columns) catalog[[column]] <- NA_character_
if (nrow(ncbi_rows)) {
  for (i in seq_len(nrow(catalog))) {
    match_index <- match(catalog$ncbi_gene_id[[i]], ncbi_rows$ncbi_gene_id)
    if (is.na(match_index)) match_index <- match(catalog$gene_symbol[[i]], ncbi_rows$ncbi_official_symbol)
    if (!is.na(match_index)) for (column in ncbi_columns) catalog[[column]][[i]] <- as.character(ncbi_rows[[column]][[match_index]])
  }
}
catalog$retrieved_at <- format(Sys.time(), tz="UTC", usetz=TRUE)

mapping <- catalog[, c("source_row","original_id","original_id_type","original_cell","token_index","input_id",
  "uniprot_accession","candidate_count","mapping_status","preferred_candidate","preferred_reason","organism_match")]
unmapped <- catalog[catalog$mapping_status == "unmapped", , drop=FALSE]
ambiguous <- catalog[catalog$mapping_status == "ambiguous", , drop=FALSE]
readr::write_csv(mapping, file.path(tables_dir, "id_mapping.csv"), na="")
readr::write_csv(catalog, file.path(tables_dir, "protein_catalog.csv"), na="")
readr::write_csv(unmapped, file.path(tables_dir, "unmapped.csv"), na="")
readr::write_csv(ambiguous, file.path(tables_dir, "ambiguous.csv"), na="")
workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "Mapping"); openxlsx::writeData(workbook, "Mapping", mapping)
openxlsx::addWorksheet(workbook, "Protein catalog"); openxlsx::writeData(workbook, "Protein catalog", catalog)
openxlsx::addWorksheet(workbook, "Unmapped"); openxlsx::writeData(workbook, "Unmapped", unmapped)
openxlsx::addWorksheet(workbook, "Ambiguous"); openxlsx::writeData(workbook, "Ambiguous", ambiguous)
openxlsx::saveWorkbook(workbook, file.path(tables_dir, "protein_mapping.xlsx"), overwrite=TRUE)

metadata <- list(run_id=parameters$`run-id`, started_at=started, finished_at=format(Sys.time(), tz="UTC", usetz=TRUE),
  organism=parameters$organism, tax_id=parameters$`tax-id`, input_identifier_type=parameters$`id-type`,
  input_identifier_column=parameters$`id-column`, uniprot_endpoint=uniprot_endpoint,
  uniprot_release=uniprot_release, ncbi_endpoint=NCBI_BASE, ncbi_status=ncbi_status,
  input_count=nrow(expanded), unique_id_count=length(unique_ids),
  expanded_group_count=sum(table(expanded$source_row) > 1),
  mapped_unique_count=length(unique(catalog$input_id[catalog$mapping_status=="mapped_unique"])),
  unmapped_count=length(unique(catalog$input_id[catalog$mapping_status=="unmapped"])),
  ambiguous_count=length(unique(catalog$input_id[catalog$mapping_status=="ambiguous"])),
  organism_mismatch_count=sum(!catalog$organism_match & catalog$mapping_status != "unmapped", na.rm=TRUE), refresh=refresh)
write_json(metadata, file.path(raw_dir, "metadata.json"))
write_json(metadata, file.path(output_root, "latest_metadata.json"))
cat(jsonlite::toJSON(metadata, auto_unbox=TRUE), "\n")
