cache_connect <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- DBI::dbConnect(RSQLite::SQLite(), path)
  DBI::dbExecute(connection, paste(
    "CREATE TABLE IF NOT EXISTS annotation_cache (",
    "id INTEGER PRIMARY KEY AUTOINCREMENT, source TEXT NOT NULL, query_id TEXT NOT NULL,",
    "query_type TEXT NOT NULL, tax_id TEXT NOT NULL, payload TEXT NOT NULL, retrieved_at TEXT NOT NULL)"))
  DBI::dbExecute(connection, "CREATE INDEX IF NOT EXISTS cache_lookup ON annotation_cache(source, query_id, query_type, tax_id, retrieved_at)")
  connection
}

cache_get <- function(connection, source, query_id, query_type, tax_id) {
  query <- DBI::dbGetQuery(connection, paste(
    "SELECT payload FROM annotation_cache WHERE source=? AND query_id=? AND query_type=? AND tax_id=?",
    "ORDER BY id DESC LIMIT 1"), params = list(source, query_id, query_type, tax_id))
  if (!nrow(query)) return(NULL)
  jsonlite::fromJSON(query$payload[[1]], simplifyVector = FALSE)
}

cache_put <- function(connection, source, query_id, query_type, tax_id, payload) {
  DBI::dbExecute(connection, paste(
    "INSERT INTO annotation_cache(source, query_id, query_type, tax_id, payload, retrieved_at)",
    "VALUES(?, ?, ?, ?, ?, ?)"), params = list(source, query_id, query_type, tax_id,
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), format(Sys.time(), tz = "UTC", usetz = TRUE)))
}

