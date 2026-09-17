`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

required_packages <- c("httr2", "jsonlite", "readr", "DBI", "RSQLite", "openxlsx")

check_dependencies <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Pacotes R ausentes: ", paste(missing, collapse = ", "))
}

parse_arguments <- function(args) {
  if (length(args) %% 2 != 0) stop("Argumentos devem usar pares --nome valor.")
  result <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    if (key == args[[i]] || !nzchar(key)) stop("Argumento inválido: ", args[[i]])
    result[[key]] <- args[[i + 1]]
  }
  required <- c("input", "output", "id-column", "id-type", "tax-id", "organism", "cache", "refresh", "run-id", "project")
  absent <- required[!required %in% names(result)]
  if (length(absent)) stop("Argumentos ausentes: ", paste(absent, collapse = ", "))
  result
}

expand_identifiers <- function(input, id_column, id_type) {
  if (!id_column %in% names(input)) stop("Coluna identificadora não encontrada: ", id_column)
  rows <- list()
  for (i in seq_len(nrow(input))) {
    original <- input[[id_column]][[i]]
    if (is.na(original) || !nzchar(trimws(as.character(original)))) next
    original <- as.character(original)
    tokens <- unlist(strsplit(trimws(original), "[;[:space:]]+"))
    tokens <- tokens[nzchar(tokens)]
    for (j in seq_along(tokens)) rows[[length(rows) + 1]] <- data.frame(
      source_row = i, original_id = original, original_id_type = id_type,
      original_cell = original, token_index = j, input_id = tokens[[j]], stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame(source_row=integer(), original_id=character(), original_id_type=character(),
    original_cell=character(), token_index=integer(), input_id=character()))
  do.call(rbind, rows)
}

collapse_values <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")

write_json <- function(value, path) jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, na = "null")

