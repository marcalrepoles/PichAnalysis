NCBI_BASE <- "https://api.ncbi.nlm.nih.gov/datasets/v2"

parse_ncbi_gene <- function(gene) {
  summaries <- gene$summary %||% list()
  summary <- if (length(summaries)) summaries[[1]] else list()
  authority <- gene$nomenclatureAuthority %||% list()
  data.frame(ncbi_gene_id=as.character(gene$geneId %||% NA_character_),
    ncbi_official_symbol=gene$symbol %||% NA_character_,
    ncbi_gene_description=gene$description %||% NA_character_, ncbi_gene_type=gene$type %||% NA_character_,
    synonyms=collapse_values(unlist(gene$synonyms %||% list())),
    hgnc_identifier=authority$identifier %||% NA_character_,
    ncbi_ensembl_gene_ids=collapse_values(unlist(gene$ensemblGeneIds %||% list())),
    swissprot_accessions=collapse_values(unlist(gene$swissProtAccessions %||% list())),
    ncbi_summary=summary$description %||% NA_character_, ncbi_summary_source=summary$source %||% NA_character_,
    ncbi_summary_date=summary$date %||% NA_character_, ncbi_tax_id=as.character(gene$taxId %||% NA_character_),
    ncbi_taxonomy_name=gene$taxname %||% NA_character_, stringsAsFactors=FALSE)
}

parse_ncbi_response <- function(payload) {
  reports <- payload$reports %||% list()
  rows <- lapply(reports, function(item) parse_ncbi_gene(item$gene %||% item))
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

fetch_ncbi_batch <- function(values, mode, tax_id) {
  encoded <- vapply(values, URLencode, character(1), reserved=TRUE)
  if (mode == "id") endpoint <- paste0(NCBI_BASE, "/gene/id/", paste(encoded, collapse=","), "/dataset_report")
  else if (mode == "accession") endpoint <- paste0(NCBI_BASE, "/gene/accession/", paste(encoded, collapse=","), "/dataset_report")
  else endpoint <- paste0(NCBI_BASE, "/gene/symbol/", paste(encoded, collapse=","), "/taxon/", tax_id, "/dataset_report")
  request <- httr2::request(endpoint) |> httr2::req_headers(Accept="application/json")
  api_key <- Sys.getenv("NCBI_API_KEY", unset="")
  if (nzchar(api_key)) request <- request |> httr2::req_url_query(api_key=api_key)
  response <- perform_request(request)
  list(data=httr2::resp_body_json(response, simplifyVector=FALSE), endpoint=endpoint)
}

