safe_name <- function(value) {
  clean <- gsub("[^A-Za-z0-9._-]+", "_", iconv(value, to="ASCII//TRANSLIT"))
  clean <- gsub("^_+|_+$", "", clean)
  ifelse(nzchar(clean), clean, "condition")
}

presence_from_values <- function(values, threshold=0, zero_is_missing=TRUE) {
  numeric <- suppressWarnings(as.numeric(values))
  present <- !is.na(numeric) & numeric > threshold
  if (!zero_is_missing && threshold == 0) present <- !is.na(numeric) & numeric >= 0
  as.integer(present)
}

reproducible <- function(n_detected, n_total, mode="count", value=2) {
  if (mode == "count") n_detected >= value else (n_detected / n_total) >= value
}

classify_two <- function(a_detected, b_detected, a_repro, b_repro, a_fraction, b_fraction, predominant=TRUE,
                         a_name="A", b_name="B") {
  strict <- ifelse(a_repro & b_detected == 0, paste0(a_name, "-specific"),
    ifelse(b_repro & a_detected == 0, paste0(b_name, "-specific"),
    ifelse(a_repro & b_repro, "Shared",
    ifelse(a_detected + b_detected > 0 & !a_repro & !b_repro, "Sporadic",
    ifelse(a_repro, paste0("Reproducible in ", a_name, " (not specific)"),
    ifelse(b_repro, paste0("Reproducible in ", b_name, " (not specific)"), "Not reproducibly detected"))))))
  exploratory <- rep(NA_character_, length(strict))
  if (predominant) {
    exploratory[a_fraction > b_fraction & b_detected > 0] <- paste0("Predominantly ", a_name)
    exploratory[b_fraction > a_fraction & a_detected > 0] <- paste0("Predominantly ", b_name)
  }
  data.frame(classification=strict, predominant_class=exploratory, stringsAsFactors=FALSE)
}

analyze_presence <- function(input, identifier_column, column_map, threshold, zero_is_missing,
                             rule_mode, rule_value, predominant=TRUE) {
  if (!identifier_column %in% names(input)) stop("Coluna identificadora ausente.")
  quantitative <- vapply(column_map, function(x) x$column, character(1))
  absent <- setdiff(quantitative, names(input))
  if (length(absent)) stop("Colunas quantitativas ausentes: ", paste(absent, collapse=", "))
  identity <- data.frame(source_row=seq_len(nrow(input)), original_id=as.character(input[[identifier_column]]),
    identifier_type=NA_character_, protein_group=grepl("[;[:space:]]", as.character(input[[identifier_column]])), stringsAsFactors=FALSE)
  binary <- as.data.frame(lapply(input[quantitative], presence_from_values,
    threshold=threshold, zero_is_missing=zero_is_missing), check.names=FALSE)
  values <- cbind(identity, input[quantitative])
  presence <- cbind(identity, binary)
  conditions <- unique(vapply(column_map, function(x) x$condition, character(1)))
  detection_rows <- list()
  for (condition in conditions) {
    columns <- vapply(Filter(function(x) identical(x$condition, condition), column_map), function(x) x$column, character(1))
    detected <- rowSums(binary[columns], na.rm=TRUE)
    total <- length(columns)
    detection_rows[[length(detection_rows)+1]] <- data.frame(identity,
      condition=condition, n_detected=detected, n_total_replicates=total,
      detection_fraction=detected/total,
      reproducible=reproducible(detected,total,rule_mode,rule_value), stringsAsFactors=FALSE)
  }
  detection <- do.call(rbind, detection_rows)
  wide_repro <- sapply(conditions, function(condition) detection$reproducible[detection$condition==condition])
  if (is.null(dim(wide_repro))) wide_repro <- matrix(wide_repro, nrow=nrow(input), ncol=length(conditions))
  pattern <- apply(wide_repro, 1, function(row) paste(paste0(conditions,"=",as.integer(row)),collapse=";"))
  classification <- identity
  classification$detection_pattern <- pattern
  classification$classification <- ifelse(rowSums(wide_repro)>0, "Reproducibly detected", "Not reproducibly detected")
  classification$predominant_class <- NA_character_
  if (length(conditions)==2) {
    a <- detection[detection$condition==conditions[[1]],]
    b <- detection[detection$condition==conditions[[2]],]
    classes <- classify_two(a$n_detected,b$n_detected,a$reproducible,b$reproducible,
      a$detection_fraction,b$detection_fraction,predominant,conditions[[1]],conditions[[2]])
    classification$classification <- classes$classification
    classification$predominant_class <- classes$predominant_class
  }
  list(identity=identity, values=values, presence=presence, detection=detection,
       classification=classification, conditions=conditions, column_map=column_map)
}
