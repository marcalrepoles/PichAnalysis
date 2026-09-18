split_mito_ids <- function(value) {
  if (!length(value) || is.na(value) || !nzchar(trimws(as.character(value)))) return(character())
  x <- unlist(strsplit(trimws(as.character(value)), "[;,|[:space:]]+"))
  unique(x[nzchar(x) & toupper(x) != "NA"])
}

validate_mitocarta_organism <- function(organism, tax_id) {
  if (!identical(as.character(organism), "Homo sapiens") || !identical(as.character(tax_id), "9606"))
    stop("Unsupported organism: MitoCarta analysis currently supports only Homo sapiens (NCBI Taxonomy ID 9606).")
  TRUE
}

empty_mito_mapping <- function() data.frame(source_row=integer(), original_id=character(), canonical_gene_key=character(),
  gene_symbol=character(), ncbi_gene_id=character(), uniprot_accession=character(), protein_name=character(),
  mapping_route=character(), mapping_status=character(), is_mitocarta=logical(), candidate_gene_symbol=character(), stringsAsFactors=FALSE)

canonical_gene_key <- function(row) {
  ncbi <- trimws(as.character(row$ncbi_gene_id %||% "")); symbol <- trimws(as.character(row$gene_symbol %||% ""))
  if (!is.na(ncbi) && nzchar(ncbi)) paste0("NCBI:", ncbi) else paste0("SYMBOL:", toupper(symbol))
}

canonical_mitocarta_mapping <- function(catalog, genes) {
  if (length(setdiff(c("source_row", "original_id"), names(catalog)))) stop("Catalog is missing required MitoCarta mapping columns.")
  for (name in c("ncbi_gene_id", "gene_symbol", "uniprot_accession", "protein_name")) if (!name %in% names(catalog)) catalog[[name]] <- NA_character_
  for (name in c("ncbi_gene_id", "gene_symbol", "uniprot_accession", "is_mitocarta")) if (!name %in% names(genes)) stop("MitoCarta genes table is missing required columns.")
  genes$gene_symbol <- toupper(trimws(as.character(genes$gene_symbol)))
  genes$ncbi_gene_id <- trimws(as.character(genes$ncbi_gene_id)); genes$uniprot_accession <- trimws(as.character(genes$uniprot_accession))
  genes$is_mitocarta <- tolower(as.character(genes$is_mitocarta)) %in% c("true", "t", "1")
  rows <- list()
  for (source in unique(catalog$source_row)) {
    group <- catalog[catalog$source_row == source, , drop=FALSE]; first <- group[1, , drop=FALSE]
    collect <- function(column) unique(unlist(lapply(group[[column]], split_mito_ids)))
    collapse <- function(column) paste(unique(as.character(group[[column]])[!is.na(group[[column]]) & nzchar(as.character(group[[column]]))]), collapse=";")
    ncbi <- collect("ncbi_gene_id"); symbols <- toupper(collect("gene_symbol")); uniprot <- collect("uniprot_accession")
    hits <- integer(); route <- "none"
    if (length(ncbi)) { hits <- which(genes$ncbi_gene_id %in% ncbi); if (length(hits)) route <- "NCBI Gene ID" }
    if (!length(hits) && length(symbols)) {
      symbol_hits <- which(genes$gene_symbol %in% symbols)
      hits <- symbol_hits
      if (length(hits)) route <- "Gene symbol"
    }
    if (!length(hits) && length(uniprot)) {
      hits <- unique(unlist(lapply(uniprot, function(up) which(vapply(genes$uniprot_accession, function(value) up %in% split_mito_ids(value), logical(1))))))
      if (length(hits)) route <- "UniProt"
    }
    hits <- unique(hits)
    status <- if (!length(hits)) "unmapped" else if (length(hits) == 1) "mapped_unique" else "ambiguous"
    if (!length(hits)) hits <- NA_integer_
    for (hit in hits) {
      gene <- if (is.na(hit)) NULL else genes[hit, , drop=FALSE]
      rows[[length(rows)+1]] <- data.frame(source_row=as.integer(first$source_row), original_id=as.character(first$original_id),
        canonical_gene_key=if (is.null(gene)) NA_character_ else canonical_gene_key(gene),
        gene_symbol=if (is.null(gene)) collapse("gene_symbol") else gene$gene_symbol,
        ncbi_gene_id=if (is.null(gene)) collapse("ncbi_gene_id") else gene$ncbi_gene_id,
        uniprot_accession=if (is.null(gene)) collapse("uniprot_accession") else gene$uniprot_accession,
        protein_name=collapse("protein_name"), mapping_route=route, mapping_status=status,
        is_mitocarta=if (is.null(gene)) NA else gene$is_mitocarta, candidate_gene_symbol=if (is.null(gene)) NA_character_ else gene$gene_symbol,
        stringsAsFactors=FALSE)
    }
  }
  if (!length(rows)) return(empty_mito_mapping())
  unique(do.call(rbind, rows))
}

build_mitocarta_membership <- function(mapping, genes) {
  selected <- mapping[mapping$mapping_status == "mapped_unique" & !is.na(mapping$canonical_gene_key), , drop=FALSE]
  columns <- c("canonical_gene_key", "gene_symbol", "ncbi_gene_id", "uniprot_accession", "is_mitocarta", "source_rows")
  if (!nrow(selected)) return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))), columns))
  keys <- unique(selected$canonical_gene_key)
  result <- do.call(rbind, lapply(keys, function(key) {
    x <- selected[selected$canonical_gene_key == key, , drop=FALSE]
    data.frame(canonical_gene_key=key, gene_symbol=x$gene_symbol[[1]], ncbi_gene_id=x$ncbi_gene_id[[1]],
      uniprot_accession=x$uniprot_accession[[1]], is_mitocarta=isTRUE(x$is_mitocarta[[1]]),
      source_rows=paste(sort(unique(x$source_row)), collapse=";"), stringsAsFactors=FALSE)
  }))
  annotation <- genes[!duplicated(toupper(genes$gene_symbol)), , drop=FALSE]
  extra <- setdiff(names(annotation), c("gene_symbol", "ncbi_gene_id", "uniprot_accession", "is_mitocarta"))
  if (length(extra)) result <- merge(result, annotation[, c("gene_symbol", extra), drop=FALSE], by="gene_symbol", all.x=TRUE, sort=FALSE)
  result
}

adjust_mito_target_background <- function(target, background, authorized=FALSE) {
  target <- unique(target[!is.na(target) & nzchar(target)]); background <- unique(background[!is.na(background) & nzchar(background)])
  outside <- setdiff(target, background)
  if (length(outside) && !authorized) return(list(ok=FALSE,target=target,background=background,outside=outside))
  list(ok=TRUE,target=intersect(target,background),background=background,outside=outside)
}

mitocarta_coverage <- function(membership, target, background) {
  t <- membership[membership$canonical_gene_key %in% target, , drop=FALSE]; b <- membership[membership$canonical_gene_key %in% background, , drop=FALSE]
  tm <- sum(t$is_mitocarta); bm <- sum(b$is_mitocarta)
  data.frame(metric=c("Target genes","Target MitoCarta genes","Target non-MitoCarta genes","Target MitoCarta fraction",
    "Background genes","Background MitoCarta genes","Background non-MitoCarta genes","Background MitoCarta fraction"),
    value=as.character(c(nrow(t),tm,nrow(t)-tm,if(nrow(t))tm/nrow(t) else 0,nrow(b),bm,nrow(b)-bm,if(nrow(b))bm/nrow(b) else 0)),stringsAsFactors=FALSE)
}

overall_mitocarta_enrichment <- function(membership, target, background, alternative="greater") {
  target_rows <- membership[membership$canonical_gene_key %in% target, , drop=FALSE]
  reference <- setdiff(background, target); reference_rows <- membership[membership$canonical_gene_key %in% reference, , drop=FALSE]
  a <- sum(target_rows$is_mitocarta); b <- nrow(target_rows)-a; c <- sum(reference_rows$is_mitocarta); d <- nrow(reference_rows)-c
  matrix_value <- matrix(c(a,c,b,d), nrow=2, dimnames=list(c("Target","Background outside target"),c("MitoCarta","non-MitoCarta")))
  fit <- tryCatch(stats::fisher.test(matrix_value, alternative=alternative), error=function(e) NULL)
  data.frame(target_mitocarta=a,target_non_mitocarta=b,reference_mitocarta=c,reference_non_mitocarta=d,
    odds_ratio=if(is.null(fit))NA_real_ else unname(fit$estimate),ci_low=if(is.null(fit))NA_real_ else fit$conf.int[[1]],
    ci_high=if(is.null(fit))NA_real_ else fit$conf.int[[2]],p_value=if(is.null(fit))NA_real_ else fit$p.value,
    alternative=alternative,stringsAsFactors=FALSE)
}

annotation_frequency <- function(annotation, target_mito, label) {
  pairs <- unique(annotation[annotation$gene_symbol %in% target_mito$gene_symbol, c(label,"gene_symbol"), drop=FALSE])
  columns <- c(label,"Gene_count","Fraction_of_MitoCarta_target","Fraction_of_all_target","Genes")
  if (!nrow(pairs)) return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))),columns))
  groups <- split(pairs$gene_symbol,pairs[[label]])
  out <- do.call(rbind,lapply(names(groups),function(name)data.frame(name=name,Gene_count=length(unique(groups[[name]])),
    Fraction_of_MitoCarta_target=if(nrow(target_mito))length(unique(groups[[name]]))/nrow(target_mito) else 0,
    Fraction_of_all_target=NA_real_,Genes=paste(sort(unique(groups[[name]])),collapse=";"),stringsAsFactors=FALSE)))
  names(out)[1] <- label;out[order(-out$Gene_count,out[[label]]),columns,drop=FALSE]
}

annotation_enrichment <- function(annotation, target_mito, background_mito, label, minimum_overlap=1L) {
  pairs <- unique(annotation[annotation$gene_symbol %in% background_mito$gene_symbol,c(label,"gene_symbol"),drop=FALSE]);terms <- unique(pairs[[label]])
  template <- data.frame(term=character(),Target_count=integer(),Target_MitoCarta_size=integer(),Background_count=integer(),Background_MitoCarta_size=integer(),GeneRatio=numeric(),BgRatio=numeric(),p_value=numeric(),FDR=numeric(),Genes=character(),stringsAsFactors=FALSE);names(template)[1]<-label
  if (!length(terms)) return(list(tested=template,excluded=template))
  rows <- do.call(rbind,lapply(terms,function(term){genes<-unique(pairs$gene_symbol[pairs[[label]]==term]);hits<-intersect(target_mito$gene_symbol,genes);data.frame(term=term,Target_count=length(hits),Target_MitoCarta_size=nrow(target_mito),Background_count=length(genes),Background_MitoCarta_size=nrow(background_mito),GeneRatio=if(nrow(target_mito))length(hits)/nrow(target_mito) else 0,BgRatio=if(nrow(background_mito))length(genes)/nrow(background_mito) else 0,p_value=NA_real_,FDR=NA_real_,Genes=paste(sort(hits),collapse=";"),stringsAsFactors=FALSE)}));names(rows)[1]<-label
  eligible <- rows$Target_count>=as.integer(minimum_overlap);tested<-rows[eligible,,drop=FALSE];excluded<-rows[!eligible,,drop=FALSE]
  if(nrow(tested)){tested$p_value<-vapply(seq_len(nrow(tested)),function(i)stats::phyper(tested$Target_count[i]-1,tested$Background_count[i],tested$Background_MitoCarta_size[i]-tested$Background_count[i],tested$Target_MitoCarta_size[i],lower.tail=FALSE),numeric(1));tested$FDR<-stats::p.adjust(tested$p_value,"BH");tested<-tested[order(tested$FDR,tested$p_value,-tested$Target_count),]}
  list(tested=tested,excluded=excluded)
}

add_hierarchy <- function(frame,hierarchy){if(!nrow(frame)||!nrow(hierarchy))return(frame);merge(frame,hierarchy,by="MitoPathway",all.x=TRUE,sort=FALSE)}

`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
