EXPERIMENTAL_EVIDENCE_CODES <- c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP")

filter_evidence <- function(annotations, mode="all") {
  if (mode == "exclude_iea") annotations[annotations$evidence_code != "IEA",,drop=FALSE]
  else if (mode == "experimental") annotations[annotations$evidence_code %in% EXPERIMENTAL_EVIDENCE_CODES,,drop=FALSE]
  else annotations
}

deduplicate_annotations <- function(annotations) {
  annotations[!duplicated(annotations[,c("entity_id","GO_ID","ontology","evidence_code")]),,drop=FALSE]
}

go_frequency <- function(annotations, total_entities) {
  unique_pairs <- unique(annotations[,c("entity_id","GO_ID","GO_term","ontology")])
  if (!nrow(unique_pairs)) return(data.frame(GO_ID=character(),Description=character(),Ontology=character(),Protein_count=integer(),Protein_fraction=numeric()))
  counts <- aggregate(entity_id~GO_ID+GO_term+ontology,unique_pairs,function(x)length(unique(x)))
  names(counts)<-c("GO_ID","Description","Ontology","Protein_count")
  counts$Protein_fraction <- counts$Protein_count/total_entities
  counts[order(-counts$Protein_count,counts$GO_ID),]
}

go_enrichment <- function(target_annotations, background_annotations, target_entities, background_entities) {
  target_entities <- intersect(unique(target_entities),unique(background_entities))
  background_entities <- unique(background_entities)
  target_pairs <- unique(target_annotations[,c("entity_id","GO_ID","GO_term","ontology")])
  background_pairs <- unique(background_annotations[,c("entity_id","GO_ID","GO_term","ontology")])
  terms <- unique(background_pairs[,c("GO_ID","GO_term","ontology")])
  rows <- lapply(seq_len(nrow(terms)),function(i) {
    term <- terms[i,]; bg_genes <- unique(background_pairs$entity_id[background_pairs$GO_ID==term$GO_ID])
    target_genes <- intersect(unique(target_pairs$entity_id[target_pairs$GO_ID==term$GO_ID]),target_entities)
    k<-length(target_genes); M<-length(bg_genes); n<-length(target_entities); N<-length(background_entities)
    data.frame(GO_ID=term$GO_ID,Description=term$GO_term,Ontology=term$ontology,
      GeneRatio=paste0(k,"/",n),BgRatio=paste0(M,"/",N),Count=k,
      pvalue=if(n&&N&&M) phyper(k-1,M,N-M,n,lower.tail=FALSE) else NA_real_,
      genes=paste(sort(target_genes),collapse=";"),stringsAsFactors=FALSE)
  })
  result<-if(length(rows))do.call(rbind,rows)else data.frame()
  if(nrow(result)){result$p.adjust<-p.adjust(result$pvalue,method="BH");result$method<-"Hypergeometric test; Benjamini-Hochberg adjustment";result<-result[order(result$p.adjust,result$pvalue),]}
  result
}

simplify_terms <- function(enrichment, annotations, cutoff=.7) {
  if(nrow(enrichment)<2)return(enrichment)
  keep<-rep(TRUE,nrow(enrichment)); gene_sets<-lapply(enrichment$GO_ID,function(go)unique(annotations$entity_id[annotations$GO_ID==go]))
  for(i in seq_len(nrow(enrichment)))if(keep[[i]]&&i<nrow(enrichment))for(j in (i+1):nrow(enrichment))if(keep[[j]]){
    union_size<-length(union(gene_sets[[i]],gene_sets[[j]])); similarity<-if(union_size)length(intersect(gene_sets[[i]],gene_sets[[j]]))/union_size else 0
    if(similarity>=cutoff)keep[[j]]<-FALSE
  }
  result<-enrichment[keep,,drop=FALSE];result$simplification_method<-"Jaccard similarity of associated experiment genes";result$simplification_cutoff<-cutoff;result$representative_rule<-"lowest BH-adjusted p-value";result
}

row_query <- function(row) {
  if (!is.na(row[["ncbi_gene_id"]]) && nzchar(as.character(row[["ncbi_gene_id"]]))) return(c("ENTREZID",strsplit(as.character(row[["ncbi_gene_id"]]),";",fixed=TRUE)[[1]][1]))
  if (!is.na(row[["gene_symbol"]]) && nzchar(as.character(row[["gene_symbol"]]))) return(c("SYMBOL",as.character(row[["gene_symbol"]])))
  if (!is.na(row[["uniprot_accession"]]) && nzchar(as.character(row[["uniprot_accession"]]))) return(c("UNIPROT",as.character(row[["uniprot_accession"]])))
  c(NA_character_,NA_character_)
}

catalog_entity_ids <- function(catalog) {
  vapply(seq_len(nrow(catalog)),function(i){q<-row_query(catalog[i,,drop=FALSE]);paste(q[[1]],q[[2]],sep=":")},character(1))
}

annotate_catalog <- function(catalog, orgdb) {
  rows<-list()
  for(i in seq_len(nrow(catalog))){
    query<-row_query(catalog[i,,drop=FALSE]);if(is.na(query[[1]]))next
    selected<-suppressMessages(AnnotationDbi::select(orgdb,keys=query[[2]],keytype=query[[1]],columns=c("ENTREZID","SYMBOL","UNIPROT","GO","EVIDENCE","ONTOLOGY")))
    selected<-selected[!is.na(selected$GO),,drop=FALSE];if(!nrow(selected))next
    entity<-ifelse(!is.na(selected$ENTREZID),paste0("ENTREZID:",selected$ENTREZID),paste(query,collapse=":"))
    rows[[length(rows)+1]]<-data.frame(source_row=catalog$source_row[[i]],input_id=catalog$original_id[[i]],
      uniprot_accession=catalog$uniprot_accession[[i]],gene_symbol=catalog$gene_symbol[[i]],
      ncbi_gene_id=catalog$ncbi_gene_id[[i]],entity_id=entity,GO_ID=selected$GO,
      ontology=selected$ONTOLOGY,evidence_code=selected$EVIDENCE,source="OrgDb direct annotation",propagation="direct",stringsAsFactors=FALSE)
  }
  if(!length(rows))return(data.frame())
  annotations<-do.call(rbind,rows)
  terms<-AnnotationDbi::select(GO.db::GO.db,keys=unique(annotations$GO_ID),keytype="GOID",columns=c("TERM","ONTOLOGY"))
  descriptions<-setNames(terms$TERM,terms$GOID);annotations$GO_term<-unname(descriptions[annotations$GO_ID])
  deduplicate_annotations(annotations)
}
