args <- commandArgs(trailingOnly=TRUE)
argument <- function(flag) {i<-match(flag,args);if(is.na(i)||i==length(args))stop(paste("Missing",flag));args[i+1]}
run <- argument("--run");provenance <- argument("--provenance")
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script <- if(length(script_arg))sub("^--file=","",script_arg[1]) else "r_scripts/10_mtdna_analysis.R"
source(file.path(dirname(script),"lib","mtdna_analysis.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required for mtDNA Evidence analysis.")
if(!requireNamespace("openxlsx",quietly=TRUE))stop("openxlsx is required for mtDNA Evidence workbook output.")
metadata <- jsonlite::fromJSON(file.path(run,"metadata.json"),simplifyVector=FALSE)
read_input <- function(name)utils::read.csv(file.path(run,"inputs",paste0(name,".csv")),stringsAsFactors=FALSE,
  check.names=FALSE,colClasses="character")
read_mapping <- function(name)utils::read.csv(file.path(run,"mapping",name),stringsAsFactors=FALSE,
  check.names=FALSE,colClasses="character")
write_table <- function(frame,relative) {
  path<-file.path(run,relative);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(frame,path,row.names=FALSE,na="")
}
target <- sort(unique(read_input("target")$entity_key))
background <- sort(unique(read_input("background")$entity_key))
if(!all(target %in% background))stop("Target outside experimental background")
entities <- read_input("entities");evidence <- read_input("evidence_records")
memberships <- read_input("category_membership");categories <- read_input("categories")
mapping <- read_mapping("experimental_entity_mapping.csv")
ambiguous <- read_mapping("ambiguous.csv");unmapped <- read_mapping("unmapped.csv")
frequency <- mtdna_frequency(categories,memberships,target,background)
sources <- mtdna_sources(evidence,target)
enrichment <- mtdna_enrichment(categories,memberships,target,background,
  as.integer(metadata$minimum_overlap),as.numeric(metadata$fdr_cutoff))
target_evidence <- evidence[evidence$entity_key %in% target,,drop=FALSE]
resolved <- mapping[nzchar(mapping$entity_key),,drop=FALSE]
entity_summary <- do.call(rbind,lapply(sort(unique(resolved$entity_key)),function(key) {
  rows <- resolved[resolved$entity_key==key,,drop=FALSE]
  part <- evidence[evidence$entity_key==key,,drop=FALSE]
  entity <- entities[entities$entity_key==key,,drop=FALSE]
  origins <- sort(unique(as.character(part$source)));cats <- sort(unique(as.character(part$evidence_category)))
  data.frame(entity_key=key,gene_symbol=if(nrow(entity))entity$gene_symbol[1] else rows$gene_symbol[1],
    ncbi_gene_id=if(nrow(entity))entity$ncbi_gene_id[1] else if(grepl("^NCBI:",key))sub("^NCBI:","",key) else "",
    experimental_accessions=paste(sort(unique(unlist(strsplit(paste(rows$original_uniprot,collapse=";"),";",fixed=TRUE)))),collapse=";"),
    source_rows=paste(sort(unique(rows$source_row)),collapse=";"),has_mtdna_related_evidence=nrow(part)>0,
    is_mtdna_encoded="NCBI mtDNA genome" %in% origins,source_count=length(origins),
    evidence_record_count=nrow(part),category_count=length(cats),sources=paste(origins,collapse=";"),
    categories=paste(cats,collapse=";"),stringsAsFactors=FALSE)
}))
classes_path <- file.path(run,"inputs","classifications.csv")
classes <- if(file.exists(classes_path))utils::read.csv(classes_path,stringsAsFactors=FALSE,check.names=FALSE) else
  data.frame(source_row=character(),classification=character())
if(nrow(classes)) {
  classes <- merge(classes,resolved[,c("source_row","entity_key"),drop=FALSE],by="source_row")
  classes <- unique(classes[,c("entity_key","classification"),drop=FALSE])
  reproducible <- classes[!classes$classification %in% c("Not reproducibly detected","Sporadic"),,drop=FALSE]
  reproducible$classification <- "Reproducibly detected"
  classes <- unique(rbind(classes,reproducible))
}
comparison <- mtdna_comparison(classes,memberships,background)
entity_to_evidence <- target_evidence[,c("entity_key","gene_symbol","source","evidence_category",
  "source_identifier","source_name","evidence_code","reference"),drop=FALSE]
category_to_entities <- unique(target_evidence[,c("evidence_category","entity_key","gene_symbol"),drop=FALSE])
names(category_to_entities)[1] <- "category"
source_to_entities <- unique(target_evidence[,c("source","entity_key","gene_symbol"),drop=FALSE])
summary <- data.frame(metric=c("Input rows","Unique resolved entities","Ambiguous entities","Unmapped entities",
  "Entities without mtDNA evidence","Target definition","Target entities","Background entities",
  "mtDNA Evidence snapshot","MitoCarta source snapshot","GO source snapshot","NCBI reference",
  "Target entities with mtDNA evidence","Target mtDNA-encoded entities","Target MitoCarta-supported entities",
  "Target GO-supported entities","Target multiple-source entities","Categories represented","Categories tested",
  "Significant categories","Minimum overlap","FDR method","FDR threshold"),
  value=as.character(c(nrow(mapping),nrow(entity_summary),nrow(ambiguous),nrow(unmapped),
    sum(!entity_summary$has_mtdna_related_evidence),metadata$target_definition,length(target),length(background),
    metadata$snapshot_id,metadata$mitocarta_snapshot_id,metadata$go_snapshot_id,metadata$ncbi_accession,
    sum(sources$agreement$source_count>0),sum(sources$agreement$mtdna_encoded),
    sum(sources$agreement$mitocarta_evidence),sum(sources$agreement$go_evidence),
    sum(sources$agreement$source_count>1),sum(frequency$target_entity_count>0),
    nrow(enrichment$all),nrow(enrichment$significant),metadata$minimum_overlap,
    "Benjamini-Hochberg",metadata$fdr_cutoff)),stringsAsFactors=FALSE)
tables <- list(
  "evidence/entity_evidence_summary.csv"=entity_summary,
  "evidence/target_evidence_records.csv"=target_evidence,
  "categories/category_frequency.csv"=frequency,
  "sources/source_frequency.csv"=sources$frequency,
  "sources/source_agreement_target.csv"=sources$agreement,
  "sources/source_count_distribution.csv"=sources$distribution,
  "enrichment/category_enrichment_all.csv"=enrichment$all,
  "enrichment/category_enrichment_significant.csv"=enrichment$significant,
  "enrichment/category_enrichment_excluded.csv"=enrichment$excluded,
  "comparison/target_category_matrix.csv"=comparison,
  "navigation/entity_to_evidence.csv"=entity_to_evidence,
  "navigation/category_to_entities.csv"=category_to_entities,
  "navigation/source_to_entities.csv"=source_to_entities,
  "summary.csv"=summary)
for(name in names(tables))write_table(tables[[name]],name)
mtdna_plots(run,frequency,sources,enrichment,target_evidence,comparison,as.integer(metadata$top_n))
workbook <- openxlsx::createWorkbook()
sheets <- list(Summary=summary,`Entity mapping`=mapping,Ambiguous=ambiguous,Unmapped=unmapped,
  `Entity evidence`=entity_summary,`Evidence records`=target_evidence,
  `Category frequency`=frequency,`Category enrichment`=enrichment$all,
  `Enrichment excluded`=enrichment$excluded,`Source frequency`=sources$frequency,
  `Source agreement`=sources$agreement,`Source count`=sources$distribution,
  `Target comparison`=comparison,`Entity to evidence`=entity_to_evidence,
  `Category to entities`=category_to_entities,`Source to entities`=source_to_entities)
for(name in names(sheets)){openxlsx::addWorksheet(workbook,name);openxlsx::writeData(workbook,name,sheets[[name]])}
openxlsx::saveWorkbook(workbook,file.path(run,"mtDNA_analysis.xlsx"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(provenance,"R_session_info.txt"))
cat(jsonlite::toJSON(list(run_id=metadata$run_id,tested=nrow(enrichment$all),
  significant=nrow(enrichment$significant)),auto_unbox=TRUE),"\n")
