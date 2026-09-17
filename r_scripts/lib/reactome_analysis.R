split_identifiers <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(trimws(as.character(value)))) return(character())
  tokens <- unlist(strsplit(trimws(as.character(value)), "[;,|[:space:]]+"))
  unique(tokens[nzchar(tokens) & toupper(tokens) != "NA"])
}

validate_reactome_organism <- function(organism, tax_id) {
  if (!identical(as.character(organism), "Homo sapiens") || !identical(as.character(tax_id), "9606"))
    stop("Unsupported organism: Reactome analysis currently supports only Homo sapiens (NCBI Taxonomy ID 9606).")
  TRUE
}

empty_reactome_mapping <- function() data.frame(
  source_row=integer(), original_id=character(), reactome_entity_key=character(),
  candidate_source_id=character(), uniprot_accession=character(), ncbi_gene_id=character(),
  gene_symbol=character(), protein_name=character(), mapping_route=character(),
  mapping_status=character(), stringsAsFactors=FALSE)

read_reactome_pathways <- function(path) {
  x <- utils::read.delim(path, header=FALSE, sep="\t", quote="", comment.char="",
    stringsAsFactors=FALSE, fill=TRUE, colClasses="character")
  if (ncol(x) < 3) stop("Invalid Reactome pathways file.")
  names(x)[1:3] <- c("Reactome_ID", "Pathway", "Species")
  unique(x[x$Species == "Homo sapiens" & grepl("^R-HSA-", x$Reactome_ID), c("Reactome_ID", "Pathway")])
}

read_reactome_mapping <- function(path, route=c("UniProt", "NCBI Gene fallback")) {
  route <- match.arg(route)
  x <- utils::read.delim(path, header=FALSE, sep="\t", quote="", comment.char="",
    stringsAsFactors=FALSE, fill=TRUE, colClasses="character")
  if (ncol(x) < 2) stop("Invalid Reactome identifier mapping file.")
  species <- if (ncol(x) >= 6) x[[6]] else if (ncol(x) >= 5) x[[5]] else rep("Homo sapiens", nrow(x))
  out <- data.frame(source_id=trimws(x[[1]]), Reactome_ID=trimws(x[[2]]), Species=species,
    mapping_route=route, stringsAsFactors=FALSE)
  unique(out[out$Species == "Homo sapiens" & grepl("^R-HSA-", out$Reactome_ID) & nzchar(out$source_id), ])
}

read_reactome_relations <- function(path, pathways) {
  x <- utils::read.delim(path, header=FALSE, sep="\t", quote="", comment.char="",
    stringsAsFactors=FALSE, fill=TRUE, colClasses="character")
  if (ncol(x) < 2) stop("Invalid Reactome pathway relation file.")
  rel <- unique(data.frame(parent_pathway_id=x[[1]], child_pathway_id=x[[2]], stringsAsFactors=FALSE))
  rel <- rel[grepl("^R-HSA-", rel$parent_pathway_id) & grepl("^R-HSA-", rel$child_pathway_id), ]
  names_map <- stats::setNames(pathways$Pathway, pathways$Reactome_ID)
  rel$parent_pathway_name <- unname(names_map[rel$parent_pathway_id])
  rel$child_pathway_name <- unname(names_map[rel$child_pathway_id])
  rel[, c("parent_pathway_id", "parent_pathway_name", "child_pathway_id", "child_pathway_name")]
}

canonical_mapping <- function(catalog, uniprot, ncbi) {
  required <- c("source_row", "original_id")
  if (length(setdiff(required, names(catalog)))) stop("Catalog is missing required Reactome mapping columns.")
  for (name in c("uniprot_accession", "ncbi_gene_id", "gene_symbol", "protein_name"))
    if (!name %in% names(catalog)) catalog[[name]] <- NA_character_
  rows <- list()
  known_uniprot <- unique(as.character(uniprot$source_id))
  known_ncbi <- unique(as.character(ncbi$source_id))
  for (source in unique(catalog$source_row)) {
    group <- catalog[catalog$source_row == source, , drop=FALSE]
    r <- group[1, , drop=FALSE]
    collect <- function(column) unique(unlist(lapply(group[[column]], split_identifiers)))
    collapse <- function(column) paste(unique(as.character(group[[column]])[!is.na(group[[column]]) & nzchar(as.character(group[[column]]))]), collapse=";")
    up <- intersect(collect("uniprot_accession"), known_uniprot)
    gene <- intersect(collect("ncbi_gene_id"), known_ncbi)
    candidates <- if (length(up)) up else gene
    route <- if (length(up)) "UniProt" else if (length(gene)) "NCBI Gene fallback" else "none"
    status <- if (!length(candidates)) "unmapped" else if (length(candidates) == 1) "mapped_unique" else "ambiguous"
    if (!length(candidates)) candidates <- NA_character_
    for (candidate in candidates) {
      key <- if (status == "unmapped") NA_character_ else paste0(if (route == "UniProt") "UP:" else "NCBI:", candidate)
      rows[[length(rows)+1]] <- data.frame(
        source_row=as.integer(r$source_row), original_id=as.character(r$original_id),
        reactome_entity_key=key, candidate_source_id=candidate,
        uniprot_accession=collapse("uniprot_accession"), ncbi_gene_id=collapse("ncbi_gene_id"),
        gene_symbol=collapse("gene_symbol"), protein_name=collapse("protein_name"),
        mapping_route=route, mapping_status=status, stringsAsFactors=FALSE)
    }
  }
  if (!length(rows)) return(empty_reactome_mapping())
  unique(do.call(rbind, rows))
}

build_pathway_membership <- function(mapping, uniprot, ncbi, pathways) {
  selected <- mapping[mapping$mapping_status == "mapped_unique" & !is.na(mapping$reactome_entity_key), , drop=FALSE]
  columns <- c("Reactome_ID", "Pathway", "reactome_entity_key", "Gene_symbol", "UniProt", "NCBI_Gene_ID", "Original_ID", "Mapping_route")
  if (!nrow(selected)) return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))), columns))
  rows <- list()
  for (i in seq_len(nrow(selected))) {
    r <- selected[i, , drop=FALSE]
    table <- if (r$mapping_route == "UniProt") uniprot else ncbi
    hits <- unique(table[table$source_id == r$candidate_source_id, c("source_id", "Reactome_ID"), drop=FALSE])
    if (!nrow(hits)) next
    rows[[length(rows)+1]] <- data.frame(
      Reactome_ID=hits$Reactome_ID,
      Pathway=pathways$Pathway[match(hits$Reactome_ID, pathways$Reactome_ID)],
      reactome_entity_key=rep(r$reactome_entity_key, nrow(hits)),
      Gene_symbol=rep(r$gene_symbol, nrow(hits)),
      UniProt=rep(r$uniprot_accession, nrow(hits)),
      NCBI_Gene_ID=rep(r$ncbi_gene_id, nrow(hits)),
      Original_ID=rep(r$original_id, nrow(hits)),
      Mapping_route=rep(r$mapping_route, nrow(hits)), stringsAsFactors=FALSE)
  }
  if (!length(rows)) return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))), columns))
  out <- do.call(rbind, rows)
  out <- out[!is.na(out$Pathway) & grepl("^R-HSA-", out$Reactome_ID), ]
  out[!duplicated(out[, c("reactome_entity_key", "Reactome_ID")]), columns, drop=FALSE]
}

hierarchy_ancestry <- function(relations) {
  nodes <- unique(c(relations$parent_pathway_id, relations$child_pathway_id))
  parents <- split(relations$parent_pathway_id, relations$child_pathway_id)
  roots <- setdiff(nodes, relations$child_pathway_id)
  walk <- function(node, trail=character()) {
    if (node %in% trail) return(list(top=character(), cycle=TRUE, depth=NA_integer_))
    direct <- unique(parents[[node]])
    if (!length(direct)) return(list(top=node, cycle=FALSE, depth=0L))
    results <- lapply(direct, function(parent) walk(parent, c(trail, node)))
    list(top=unique(unlist(lapply(results, `[[`, "top"))),
      cycle=any(vapply(results, `[[`, logical(1), "cycle")),
      depth=if (all(vapply(results, function(x) is.na(x$depth), logical(1)))) NA_integer_ else
        1L + min(unlist(lapply(results, `[[`, "depth")), na.rm=TRUE))
  }
  rows <- lapply(nodes, function(node) {
    result <- walk(node)
    data.frame(pathway_id=node, top_level_ancestors=paste(sort(result$top), collapse=";"),
      minimum_depth=result$depth, cycle_detected=result$cycle, stringsAsFactors=FALSE)
  })
  do.call(rbind, rows)
}

select_reactome_target <- function(catalog, definition, classification=NULL, manual_rows=NULL) {
  if (definition %in% c("All mapped entities", "all_mapped")) return(catalog)
  if (definition %in% c("Manual selection", "manual")) return(catalog[catalog$source_row %in% (manual_rows %||% integer()), , drop=FALSE])
  if (is.null(classification)) stop("Presence/Absence outputs are required for the selected Reactome target.")
  label <- if (startsWith(definition, "class:")) substring(definition, 7) else definition
  if (definition == "Reproducibly detected") {
    rows <- classification$source_row[classification$classification != "Not reproducibly detected"]
  } else if (definition == "Shared") {
    rows <- classification$source_row[classification$classification == "Shared"]
  } else if (definition == "Sporadic") {
    rows <- classification$source_row[classification$classification == "Sporadic"]
  } else {
    rows <- classification$source_row[classification$classification == label]
  }
  catalog[catalog$source_row %in% rows, , drop=FALSE]
}

adjust_target_background <- function(target, background, authorized=FALSE) {
  target <- unique(target[!is.na(target) & nzchar(target)])
  background <- unique(background[!is.na(background) & nzchar(background)])
  outside <- setdiff(target, background)
  if (length(outside) && !authorized) return(list(ok=FALSE, target=target, background=background, outside=outside))
  list(ok=TRUE, target=intersect(target, background), background=background, outside=outside)
}

reactome_frequency <- function(membership, target_entities=NULL) {
  if (is.null(target_entities)) target_entities <- unique(membership$reactome_entity_key)
  target_entities <- unique(target_entities[!is.na(target_entities)])
  pairs <- unique(membership[membership$reactome_entity_key %in% target_entities,
    c("reactome_entity_key", "Reactome_ID", "Pathway"), drop=FALSE])
  columns <- c("Reactome_ID", "Pathway", "Protein_count", "Protein_fraction")
  if (!nrow(pairs)) return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))), columns))
  out <- stats::aggregate(reactome_entity_key ~ Reactome_ID + Pathway, pairs,
    function(x) length(unique(x)))
  names(out)[3] <- "Protein_count"
  out$Protein_fraction <- out$Protein_count / length(target_entities)
  out[order(-out$Protein_count, out$Reactome_ID), columns]
}

reactome_enrichment <- function(membership, target_entities=NULL, background_entities=NULL, minimum_overlap=3L) {
  if (is.null(target_entities) && is.data.frame(membership) && is.data.frame(background_entities)) stop("Invalid Reactome enrichment arguments.")
  target_entities <- unique(target_entities[!is.na(target_entities)])
  background_entities <- unique(background_entities[!is.na(background_entities)])
  pairs <- unique(membership[membership$reactome_entity_key %in% background_entities,
    c("reactome_entity_key", "Reactome_ID", "Pathway"), drop=FALSE])
  ids <- unique(pairs$Reactome_ID)
  all_rows <- lapply(ids, function(id) {
    term <- pairs[pairs$Reactome_ID == id, , drop=FALSE]
    bg <- unique(term$reactome_entity_key)
    hit <- intersect(target_entities, bg)
    data.frame(Reactome_ID=id, Pathway=term$Pathway[[1]], Target_count=length(hit),
      Target_size=length(target_entities), Background_count=length(bg), Background_size=length(background_entities),
      GeneRatio=if (length(target_entities)) length(hit)/length(target_entities) else 0,
      BgRatio=if (length(background_entities)) length(bg)/length(background_entities) else 0,
      p_value=NA_real_, FDR=NA_real_, Entities=paste(sort(hit), collapse=";"),
      stringsAsFactors=FALSE)
  })
  template <- data.frame(Reactome_ID=character(),Pathway=character(),Target_count=integer(),Target_size=integer(),
    Background_count=integer(),Background_size=integer(),GeneRatio=numeric(),BgRatio=numeric(),p_value=numeric(),
    FDR=numeric(),Entities=character(),stringsAsFactors=FALSE)
  if (!length(all_rows)) return(list(tested=template, excluded=template))
  all <- do.call(rbind, all_rows)
  eligible <- all$Target_count >= as.integer(minimum_overlap)
  tested <- all[eligible, , drop=FALSE]
  excluded <- all[!eligible, , drop=FALSE]
  if (nrow(tested)) {
    tested$p_value <- vapply(seq_len(nrow(tested)), function(i) stats::phyper(
      tested$Target_count[i]-1, tested$Background_count[i],
      tested$Background_size[i]-tested$Background_count[i], tested$Target_size[i], lower.tail=FALSE), numeric(1))
    tested$FDR <- stats::p.adjust(tested$p_value, method="BH")
    tested <- tested[order(tested$FDR, tested$p_value, -tested$Target_count), ]
  }
  list(tested=tested, excluded=excluded)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
