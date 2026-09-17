args_all <- commandArgs(trailingOnly=FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=", args_all)]), winslash="/", mustWork=TRUE)
script_root <- dirname(script_path)
source(file.path(script_root,"lib","common.R"))
source(file.path(script_root,"lib","presence_absence.R"))
needed <- c("jsonlite","readr","openxlsx","ggplot2")
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Pacotes R ausentes: ",paste(missing,collapse=", "))
parameters <- parse_arguments(commandArgs(trailingOnly=TRUE), c("input","output","project","identifier-column","identifier-type",
  "quantification-type","column-map","conditions","zero-is-missing","threshold","rule-mode","rule-value","predominant","run-id"))
started <- format(Sys.time(),tz="UTC",usetz=TRUE)
output_root <- normalizePath(parameters$output,winslash="/",mustWork=FALSE)
run_root <- file.path(output_root,"runs",parameters$`run-id`)
snapshot <- file.path(parameters$project,"scripts","runs",paste0(parameters$`run-id`,"_presence_absence"))
if (dir.exists(run_root)||dir.exists(snapshot)) stop("run_id já existe; histórico não será sobrescrito.")
run_tables <- file.path(run_root,"tables"); run_graphs <- file.path(run_root,"graphs")
latest_tables <- file.path(output_root,"tables"); latest_graphs <- file.path(output_root,"graphs")
for (path in c(run_tables,run_graphs,latest_tables,latest_graphs,file.path(output_root,"raw"),snapshot)) dir.create(path,recursive=TRUE,showWarnings=FALSE)
file.copy(script_path,file.path(snapshot,basename(script_path)),overwrite=FALSE)
dir.create(file.path(snapshot,"lib"),showWarnings=FALSE)
for (name in c("common.R","presence_absence.R")) file.copy(file.path(script_root,"lib",name),file.path(snapshot,"lib",name),overwrite=FALSE)
write_json(parameters,file.path(snapshot,"parameters.json")); capture.output(sessionInfo(),file=file.path(snapshot,"session_info.txt"))

column_map <- jsonlite::fromJSON(parameters$`column-map`,simplifyVector=FALSE)
conditions <- unlist(jsonlite::fromJSON(parameters$conditions))
input <- readr::read_csv(parameters$input,show_col_types=FALSE,name_repair="minimal")
result <- analyze_presence(input,parameters$`identifier-column`,column_map,as.numeric(parameters$threshold),
  identical(parameters$`zero-is-missing`,"true"),parameters$`rule-mode`,as.numeric(parameters$`rule-value`),
  identical(parameters$predominant,"true"))
result$identity$identifier_type <- parameters$`identifier-type`
for (frame_name in c("values","presence","detection","classification")) result[[frame_name]]$identifier_type <- parameters$`identifier-type`

catalog_path <- file.path(parameters$project,"mapping","tables","protein_catalog.csv")
if (file.exists(catalog_path)) {
  catalog <- readr::read_csv(catalog_path,show_col_types=FALSE)
  keep <- intersect(c("source_row","uniprot_accession","gene_symbol","protein_name","mapping_status","ncbi_summary","preferred_candidate"),names(catalog))
  catalog <- catalog[,keep,drop=FALSE]
  if ("preferred_candidate" %in% names(catalog)) catalog <- catalog[order(!catalog$preferred_candidate),]
  catalog <- catalog[!duplicated(catalog$source_row),setdiff(names(catalog),"preferred_candidate"),drop=FALSE]
  result$classification <- merge(result$classification,catalog,by="source_row",all.x=TRUE,sort=FALSE)
}

summary_rows <- lapply(result$conditions,function(condition) {
  rows <- result$detection[result$detection$condition==condition,]
  data.frame(condition=condition,total_replicates=unique(rows$n_total_replicates),
    proteins_detected_any=sum(rows$n_detected>0),proteins_reproducible=sum(rows$reproducible),stringsAsFactors=FALSE)
})
condition_summary <- do.call(rbind,summary_rows)
write_set <- list(presence_matrix=result$presence,quantitative_values=result$values,
  condition_detection=result$detection,classification=result$classification,condition_summary=condition_summary)
for (name in names(write_set)) readr::write_csv(write_set[[name]],file.path(run_tables,paste0(name,".csv")),na="")

specific_files <- list()
if (length(result$conditions)==2) {
  a <- result$conditions[[1]]; b <- result$conditions[[2]]
  specific_files[[paste0(safe_name(a),"_specific")]] <- result$classification[result$classification$classification==paste0(a,"-specific"),]
  specific_files[[paste0(safe_name(b),"_specific")]] <- result$classification[result$classification$classification==paste0(b,"-specific"),]
  specific_files$shared <- result$classification[result$classification$classification=="Shared",]
  specific_files$sporadic <- result$classification[result$classification$classification=="Sporadic",]
  for (name in names(specific_files)) readr::write_csv(specific_files[[name]],file.path(run_tables,paste0(name,".csv")),na="")
}

save_plot <- function(plot,name,width=9,height=6) {
  ggplot2::ggsave(file.path(run_graphs,paste0(name,".png")),plot,width=width,height=height,dpi=160)
  ggplot2::ggsave(file.path(run_graphs,paste0(name,".pdf")),plot,width=width,height=height)
}
sample_labels <- vapply(column_map,function(x) paste(x$condition,x$replicate,sep=" / "),character(1))
sample_counts <- data.frame(sample=sample_labels,condition=vapply(column_map,function(x)x$condition,character(1)),
  replicate=vapply(column_map,function(x)x$replicate,character(1)),count=colSums(result$presence[vapply(column_map,function(x)x$column,character(1))]))
p1 <- ggplot2::ggplot(sample_counts,ggplot2::aes(sample,count,fill=condition))+ggplot2::geom_col()+
  ggplot2::labs(x="Condição / réplica",y="Proteínas detectadas",title="Proteínas detectadas por réplica")+
  ggplot2::theme_minimal()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
save_plot(p1,"protein_count_per_sample")

binary <- result$presence[vapply(column_map,function(x)x$column,character(1))]
patterns <- apply(binary,1,paste0,collapse="")
pattern_counts <- as.data.frame(table(patterns),stringsAsFactors=FALSE); names(pattern_counts)<-c("pattern","count")
p2 <- ggplot2::ggplot(pattern_counts,ggplot2::aes(reorder(pattern,count),count))+
  ggplot2::geom_col(fill="#3B82F6")+ggplot2::coord_flip()+
  ggplot2::labs(x="Padrão binário de pertencimento (ordem das réplicas)",y="Interseção",title="UpSet — interseções entre réplicas")+ggplot2::theme_minimal()
save_plot(p2,"upset_intersections")

heatmap_note <- "all rows"
heat_rows <- seq_len(nrow(binary))
if (length(heat_rows)>200) { heat_rows <- order(rowSums(binary),decreasing=TRUE)[seq_len(200)]; heatmap_note <- "top 200 by detections" }
long <- data.frame(source_row=rep(result$identity$source_row[heat_rows],times=ncol(binary)),
  sample=rep(sample_labels,each=length(heat_rows)),presence=unlist(binary[heat_rows,,drop=FALSE],use.names=FALSE))
p3 <- ggplot2::ggplot(long,ggplot2::aes(sample,factor(source_row),fill=factor(presence)))+ggplot2::geom_tile()+
  ggplot2::scale_fill_manual(values=c("0"="white","1"="#2563EB"),name="Detectada")+
  ggplot2::labs(x="Réplica",y="Linha experimental",title="Heatmap binário de presença/ausência")+ggplot2::theme_minimal()+
  ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=45,hjust=1))
save_plot(p3,"presence_heatmap",10,7)

if (length(result$conditions)==2) {
  a <- result$detection[result$detection$condition==result$conditions[[1]],]
  b <- result$detection[result$detection$condition==result$conditions[[2]],]
  scatter <- data.frame(a=a$detection_fraction,b=b$detection_fraction,class=result$classification$classification)
  p4 <- ggplot2::ggplot(scatter,ggplot2::aes(a,b,color=class))+ggplot2::geom_point(alpha=.7)+
    ggplot2::labs(x=paste("Fração detectada",result$conditions[[1]]),y=paste("Fração detectada",result$conditions[[2]]),
      title="Fração de detecção entre condições",color="Classificação")+ggplot2::theme_minimal()
  save_plot(p4,"detection_fraction")
}

workbook <- openxlsx::createWorkbook()
sheet_data <- c(list(`Presence matrix`=result$presence,`Condition detection`=result$detection,
  Classification=result$classification,`Condition summary`=condition_summary),specific_files)
for (name in names(sheet_data)) { sheet <- substr(gsub("_"," ",name),1,31); openxlsx::addWorksheet(workbook,sheet); openxlsx::writeData(workbook,sheet,sheet_data[[name]]) }
openxlsx::saveWorkbook(workbook,file.path(run_tables,"presence_absence.xlsx"),overwrite=FALSE)

unlink(latest_tables,recursive=TRUE); unlink(latest_graphs,recursive=TRUE); dir.create(latest_tables,recursive=TRUE); dir.create(latest_graphs,recursive=TRUE)
file.copy(list.files(run_tables,full.names=TRUE),latest_tables,overwrite=TRUE)
file.copy(list.files(run_graphs,full.names=TRUE),latest_graphs,overwrite=TRUE)
counts <- table(result$classification$classification)
metadata <- list(run_id=parameters$`run-id`,started_at=started,finished_at=format(Sys.time(),tz="UTC",usetz=TRUE),
  quantification_type=parameters$`quantification-type`,conditions=result$conditions,zero_is_missing=identical(parameters$`zero-is-missing`,"true"),
  threshold=as.numeric(parameters$threshold),rule_mode=parameters$`rule-mode`,rule_value=as.numeric(parameters$`rule-value`),
  total_entities=nrow(input),classification_counts=as.list(counts),heatmap_selection=heatmap_note,
  input_file=parameters$input,identifier_column=parameters$`identifier-column`)
write_json(metadata,file.path(run_root,"metadata.json")); write_json(metadata,file.path(output_root,"latest_metadata.json"))
cat(jsonlite::toJSON(metadata,auto_unbox=TRUE),"\n")
