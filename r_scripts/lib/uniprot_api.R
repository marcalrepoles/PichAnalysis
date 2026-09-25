UNIPROT_BASE <- "https://rest.uniprot.org"
UNIPROT_FIELDS_ENDPOINT <- paste0(UNIPROT_BASE, "/configure/idmapping/fields")
UNIPROT_DATABASES <- c(uniprot="UniProtKB_AC-ID", gene_symbol="Gene_Name", entrez="GeneID",
  ensembl_gene="Ensembl", ensembl_protein="Ensembl_Protein", refseq_protein="RefSeq_Protein")

validate_uniprot_database <- function(id_type) {
  response <- perform_request(httr2::request(UNIPROT_FIELDS_ENDPOINT))
  payload <- httr2::resp_body_json(response, simplifyVector=FALSE)
  names_available <- unlist(lapply(payload$groups %||% list(), function(group)
    vapply(group$items %||% list(), function(item) if (isTRUE(item$from)) item$name %||% "" else "", character(1))))
  expected <- unname(UNIPROT_DATABASES[[id_type]])
  if (is.null(expected) || !expected %in% names_available) stop("Source database is not accepted by the current UniProt API: ", expected)
  payload
}

perform_request <- function(request, attempts = 5) {
  request |>
    httr2::req_timeout(60) |>
    httr2::req_retry(max_tries = attempts, backoff = function(n) min(30, 2^(n - 1)),
      is_transient = function(resp) httr2::resp_status(resp) == 429 || httr2::resp_status(resp) >= 500) |>
    httr2::req_perform()
}

submit_uniprot_job <- function(ids, id_type, tax_id) {
  from <- unname(UNIPROT_DATABASES[[id_type]])
  if (is.null(from)) stop("Identifier type is not supported by UniProt: ", id_type)
  body <- list(from = from, to = "UniProtKB", ids = paste(ids, collapse = ","))
  if (identical(id_type, "gene_symbol")) body$taxId <- tax_id
  response <- perform_request(httr2::request(paste0(UNIPROT_BASE, "/idmapping/run")) |>
    httr2::req_method("POST") |> httr2::req_body_form(!!!body))
  httr2::resp_body_json(response, simplifyVector = TRUE)$jobId
}

wait_uniprot_job <- function(job_id, max_polls = 60) {
  for (attempt in seq_len(max_polls)) {
    response <- perform_request(httr2::request(paste0(UNIPROT_BASE, "/idmapping/status/", job_id)))
    status <- httr2::resp_body_json(response, simplifyVector = TRUE)
    if (!is.null(status$results) || !is.null(status$failedIds) || identical(status$jobStatus, "FINISHED")) return(TRUE)
    if (!is.null(status$jobStatus) && !status$jobStatus %in% c("NEW", "RUNNING")) stop("UniProt job falhou: ", status$jobStatus)
    Sys.sleep(min(10, attempt))
  }
  stop("Tempo limite de polling do UniProt excedido.")
}

fetch_uniprot_results <- function(job_id) {
  details_response <- perform_request(httr2::request(paste0(UNIPROT_BASE, "/idmapping/details/", job_id)))
  details <- httr2::resp_body_json(details_response, simplifyVector = TRUE)
  redirect <- details$redirectURL
  stream <- if (grepl("/idmapping/results/", redirect, fixed=TRUE)) sub("/results/", "/stream/", redirect, fixed=TRUE) else sub("/results/", "/results/stream/", redirect, fixed=TRUE)
  separator <- if (grepl("?", stream, fixed=TRUE)) "&" else "?"
  response <- perform_request(httr2::request(paste0(stream, separator, "format=json")))
  list(data = httr2::resp_body_json(response, simplifyVector = FALSE),
       release = httr2::resp_header(response, "x-uniprot-release"), endpoint = stream)
}

safe_values <- function(items, field) {
  if (is.null(items)) return(character())
  vapply(items, function(item) as.character(item[[field]] %||% ""), character(1))
}

comment_text <- function(comments, type) {
  selected <- Filter(function(item) identical(item$commentType, type), comments %||% list())
  if (!length(selected)) return(NA_character_)
  values <- unlist(lapply(selected, function(item) {
    direct <- safe_values(item$texts, "value")
    locations <- unlist(lapply(item$subcellularLocations %||% list(), function(loc) loc$location$value %||% ""))
    c(direct, locations)
  }))
  collapse_values(as.character(values))
}

xref_values <- function(refs, database) collapse_values(vapply(Filter(
  function(item) identical(item$database, database), refs %||% list()),
  function(item) as.character(item$id %||% ""), character(1)))

parse_uniprot_entry <- function(from, entry) {
  genes <- entry$genes %||% list()
  first_gene <- if (length(genes)) genes[[1]] else list()
  protein <- entry$proteinDescription$recommendedName$fullName$value %||%
    entry$proteinDescription$submissionNames[[1]]$fullName$value %||% NA_character_
  data.frame(input_id=from, uniprot_accession=entry$primaryAccession %||% NA_character_,
    uniprot_entry_name=entry$uniProtkbId %||% NA_character_,
    uniprot_reviewed=grepl("Swiss-Prot", entry$entryType %||% "", fixed=TRUE),
    gene_symbol=first_gene$geneName$value %||% NA_character_,
    uniprot_gene_names=collapse_values(c(first_gene$geneName$value %||% "", safe_values(first_gene$synonyms, "value"))),
    organism_name=entry$organism$scientificName %||% NA_character_, tax_id=as.character(entry$organism$taxonId %||% NA_character_),
    protein_length=entry$sequence$length %||% NA_integer_, protein_name=protein,
    uniprot_function=comment_text(entry$comments, "FUNCTION"),
    uniprot_subcellular_location=comment_text(entry$comments, "SUBCELLULAR LOCATION"),
    ncbi_gene_id=xref_values(entry$uniProtKBCrossReferences, "GeneID"),
    ensembl_gene_ids=xref_values(entry$uniProtKBCrossReferences, "Ensembl"),
    refseq_ids=xref_values(entry$uniProtKBCrossReferences, "RefSeq"), stringsAsFactors=FALSE)
}

parse_uniprot_response <- function(payload) {
  rows <- lapply(payload$results %||% list(), function(item) parse_uniprot_entry(as.character(item$from), item$to))
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

resolve_candidate_flags <- function(candidates, tax_id) {
  count <- nrow(candidates)
  status <- if (count == 0) "unmapped" else if (count == 1) "mapped_unique" else "ambiguous"
  reviewed_count <- if (count) sum(candidates$uniprot_reviewed %in% TRUE) else 0
  list(candidate_count=count, mapping_status=status,
    preferred_candidate=if (count) count > 1 & reviewed_count == 1 & candidates$uniprot_reviewed %in% TRUE else logical(),
    organism_match=if (count) !is.na(candidates$tax_id) & candidates$tax_id == tax_id else logical())
}
