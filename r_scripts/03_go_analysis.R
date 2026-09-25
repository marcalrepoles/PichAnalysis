args_all<-commandArgs(trailingOnly=FALSE);script_path<-normalizePath(sub("^--file=","",args_all[grepl("^--file=",args_all)]),winslash="/",mustWork=TRUE)
script_root<-dirname(script_path);source(file.path(script_root,"lib","common.R"));source(file.path(script_root,"lib","go_analysis.R"))
required_packages<-c("AnnotationDbi","GO.db","ggplot2","openxlsx","jsonlite","readr")
missing<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))stop("Missing R packages: ",paste(missing,collapse=", "))
required_args<-c("input","catalog","output","project","target-file","background-file","organism","tax-id","ontologies","evidence-filter","fdr-cutoff","p-cutoff","min-count","top-n","simplify","simplify-cutoff","run-id","target-name","background-name","outside-count")
p<-parse_arguments(commandArgs(trailingOnly=TRUE),required_args);started<-format(Sys.time(),tz="UTC",usetz=TRUE)
orgdb_name<-if(p$`tax-id`=="9606")"org.Hs.eg.db"else if(p$`tax-id`=="10090")"org.Mm.eg.db"else stop("Local GO is not configured for tax_id ",p$`tax-id`)
if(!requireNamespace(orgdb_name,quietly=TRUE))stop("Missing R package: ",orgdb_name)
suppressPackageStartupMessages(library(orgdb_name,character.only=TRUE));orgdb<-get(orgdb_name)
ontologies<-unlist(jsonlite::fromJSON(p$ontologies));root<-normalizePath(p$output,winslash="/",mustWork=FALSE)
run_root<-file.path(root,"runs",p$`run-id`);snapshot<-file.path(p$project,"scripts","runs",paste0(p$`run-id`,"_go_analysis"))
if(dir.exists(run_root)||dir.exists(snapshot))stop("GO run_id already exists.")
run_ann<-file.path(run_root,"annotation");run_freq<-file.path(run_root,"frequency");run_enrich<-file.path(run_root,"enrichment");run_graphs<-file.path(run_root,"graphs")
for(path in c(run_ann,run_freq,run_enrich,run_graphs,snapshot))dir.create(path,recursive=TRUE,showWarnings=FALSE)
file.copy(script_path,file.path(snapshot,basename(script_path)),overwrite=FALSE);dir.create(file.path(snapshot,"lib"),showWarnings=FALSE)
for(name in c("common.R","go_analysis.R"))file.copy(file.path(script_root,"lib",name),file.path(snapshot,"lib",name),overwrite=FALSE)
write_json(p,file.path(snapshot,"parameters.json"));capture.output(sessionInfo(),file=file.path(snapshot,"session_info.txt"))
target<-readr::read_csv(p$`target-file`,show_col_types=FALSE);background<-readr::read_csv(p$`background-file`,show_col_types=FALSE)
target_all<-annotate_catalog(target,orgdb);background_all<-annotate_catalog(background,orgdb)
target_annotations<-filter_evidence(target_all,p$`evidence-filter`);background_annotations<-filter_evidence(background_all,p$`evidence-filter`)
target_annotations<-target_annotations[target_annotations$ontology%in%ontologies,,drop=FALSE];background_annotations<-background_annotations[background_annotations$ontology%in%ontologies,,drop=FALSE]
readr::write_csv(target_annotations,file.path(run_ann,"go_annotations.csv"),na="")
annotated_rows<-unique(target_annotations$source_row);unannotated<-target[!target$source_row%in%annotated_rows,,drop=FALSE]
readr::write_csv(unannotated,file.path(run_ann,"unannotated_proteins.csv"),na="")
target_entities<-unique(c(target_all$entity_id,catalog_entity_ids(target[!target$source_row%in%unique(target_all$source_row),,drop=FALSE])))
background_entities<-unique(c(background_all$entity_id,catalog_entity_ids(background[!background$source_row%in%unique(background_all$source_row),,drop=FALSE])))
target_entities<-target_entities[!grepl("^NA:",target_entities)];background_entities<-background_entities[!grepl("^NA:",background_entities)]
summaries<-list();tables_for_book<-list(Annotations=target_annotations,Unannotated=unannotated)
save_plot<-function(plot,name){ggplot2::ggsave(file.path(run_graphs,paste0(name,".png")),plot,width=10,height=7,dpi=160);ggplot2::ggsave(file.path(run_graphs,paste0(name,".pdf")),plot,width=10,height=7)}
for(ontology in ontologies){
  ann<-target_annotations[target_annotations$ontology==ontology,,drop=FALSE];bg<-background_annotations[background_annotations$ontology==ontology,,drop=FALSE]
  frequency<-go_frequency(ann,length(unique(target_entities)));readr::write_csv(frequency,file.path(run_freq,paste0("go_",tolower(ontology),"_frequency.csv")),na="")
  enrichment<-go_enrichment(ann,bg,target_entities,background_entities)
  if(nrow(enrichment)){enrichment$passes_min_count<-enrichment$Count>=as.integer(p$`min-count`);enrichment$significant<-enrichment$passes_min_count&enrichment$pvalue<=as.numeric(p$`p-cutoff`)&enrichment$p.adjust<=as.numeric(p$`fdr-cutoff`)}
  significant<-if(nrow(enrichment))enrichment[enrichment$significant,,drop=FALSE]else enrichment
  readr::write_csv(enrichment,file.path(run_enrich,paste0("go_",tolower(ontology),"_all.csv")),na="");readr::write_csv(significant,file.path(run_enrich,paste0("go_",tolower(ontology),"_significant.csv")),na="")
  simplified<-data.frame();if(identical(p$simplify,"true")&&nrow(significant)){simplified<-simplify_terms(significant,ann,as.numeric(p$`simplify-cutoff`));readr::write_csv(simplified,file.path(run_enrich,paste0("go_",tolower(ontology),"_simplified.csv")),na="")}
  summaries[[ontology]]<-data.frame(ontology=ontology,annotated_proteins=length(unique(ann$entity_id)),terms_tested=nrow(enrichment),significant_terms=nrow(significant))
  tables_for_book[[paste(ontology,"frequency")]]<-frequency;tables_for_book[[paste(ontology,"enrichment")]]<-enrichment
  topn<-as.integer(p$`top-n`)
  if(nrow(frequency)){plotdata<-head(frequency,topn);q<-ggplot2::ggplot(plotdata,ggplot2::aes(reorder(Description,Protein_count),Protein_count))+ggplot2::geom_col(fill="#2563EB")+ggplot2::coord_flip()+ggplot2::labs(x="GO term",y="Number of proteins",title=paste("GO frequency",ontology,"(not enrichment)"))+ggplot2::theme_minimal();save_plot(q,paste0("go_",tolower(ontology),"_frequency"))}
  if(nrow(significant)){plotdata<-head(significant,topn);ratio<-vapply(strsplit(plotdata$GeneRatio,"/",fixed=TRUE),function(x)as.numeric(x[1])/as.numeric(x[2]),numeric(1));plotdata$ratio<-ratio
    q<-ggplot2::ggplot(plotdata,ggplot2::aes(ratio,reorder(Description,ratio),size=Count,color=p.adjust))+ggplot2::geom_point()+ggplot2::scale_color_viridis_c(direction=-1)+ggplot2::labs(x="GeneRatio",y="GO term",title=paste("GO enrichment",ontology),color="FDR (BH)")+ggplot2::theme_minimal();save_plot(q,paste0("go_",tolower(ontology),"_enrichment_dot"))
    plotdata$score<--log10(pmax(plotdata$p.adjust,.Machine$double.xmin));q2<-ggplot2::ggplot(plotdata,ggplot2::aes(reorder(Description,score),score))+ggplot2::geom_col(fill="#7C3AED")+ggplot2::coord_flip()+ggplot2::labs(x="GO term",y="-log10(FDR BH)",title=paste("GO enrichment",ontology))+ggplot2::theme_minimal();save_plot(q2,paste0("go_",tolower(ontology),"_enrichment_bar"))}
}
summary<-do.call(rbind,summaries);summary$target_proteins<-nrow(target);summary$unannotated_proteins<-nrow(unannotated);summary$background_proteins<-nrow(background)
readr::write_csv(summary,file.path(run_root,"summary.csv"));tables_for_book<-c(list(Summary=summary),tables_for_book)
wb<-openxlsx::createWorkbook();for(name in names(tables_for_book)){sheet<-substr(name,1,31);openxlsx::addWorksheet(wb,sheet);openxlsx::writeData(wb,sheet,tables_for_book[[name]])};openxlsx::saveWorkbook(wb,file.path(run_root,"GO_analysis.xlsx"),overwrite=FALSE)
for(folder in c("annotation","frequency","enrichment","graphs")){latest<-file.path(root,folder);unlink(latest,recursive=TRUE);dir.create(latest,recursive=TRUE);file.copy(list.files(file.path(run_root,folder),full.names=TRUE),latest,overwrite=TRUE)}
file.copy(file.path(run_root,"summary.csv"),file.path(root,"summary.csv"),overwrite=TRUE);file.copy(file.path(run_root,"GO_analysis.xlsx"),file.path(root,"GO_analysis.xlsx"),overwrite=TRUE)
metadata<-list(run_id=p$`run-id`,started_at=started,finished_at=format(Sys.time(),tz="UTC",usetz=TRUE),organism=p$organism,tax_id=p$`tax-id`,orgdb=orgdb_name,orgdb_version=as.character(packageVersion(orgdb_name)),go_db_version=as.character(packageVersion("GO.db")),r_version=R.version.string,evidence_filter=p$`evidence-filter`,ontologies=ontologies,target_name=p$`target-name`,background_name=p$`background-name`,target_count=nrow(target),background_count=nrow(background),target_outside_background=as.integer(p$`outside-count`),annotated_count=length(unique(target_annotations$source_row)),unannotated_count=nrow(unannotated),enrichment_method="Hypergeometric test",p_adjust_method="Benjamini-Hochberg",fdr_cutoff=as.numeric(p$`fdr-cutoff`),min_count=as.integer(p$`min-count`),simplify=identical(p$simplify,"true"),simplify_method=if(identical(p$simplify,"true"))"Jaccard similarity of experiment genes"else NA)
write_json(metadata,file.path(run_root,"metadata.json"));write_json(metadata,file.path(root,"latest_metadata.json"));cat(jsonlite::toJSON(metadata,auto_unbox=TRUE),"\n")
