args <- commandArgs(trailingOnly=TRUE)
argument <- function(flag) {i <- match(flag,args);if(is.na(i)||i==length(args))stop(paste("Missing",flag));args[i+1]}
run <- argument("--run");provenance <- argument("--provenance")
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script <- if(length(script_arg))sub("^--file=","",script_arg[1]) else "r_scripts/09_complex_analysis.R"
source(file.path(dirname(script),"lib","complex_analysis.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("The jsonlite package is required for Complex Portal analysis.")
if(!requireNamespace("openxlsx",quietly=TRUE))stop("The openxlsx package is required for Complex Portal workbook output.")
metadata <- jsonlite::fromJSON(file.path(run,"metadata.json"),simplifyVector=FALSE)
read_input <- function(name)utils::read.csv(file.path(run,"inputs",paste0(name,".csv")),stringsAsFactors=FALSE,check.names=FALSE,colClasses="character")
write_table <- function(frame,relative){path<-file.path(run,relative);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);utils::write.csv(frame,path,row.names=FALSE,na="")}
complexes <- read_input("complexes");expanded <- read_input("expanded_protein_components")
direct <- read_input("direct_participants");nonprotein <- read_input("nonprotein_participants");nested <- read_input("nested_complexes")
target <- unique(read_input("target")$canonical_uniprot);background <- unique(read_input("background")$canonical_uniprot)
mapping <- utils::read.csv(file.path(run,"mapping/protein_mapping.csv"),stringsAsFactors=FALSE,check.names=FALSE)
ambiguous <- utils::read.csv(file.path(run,"mapping/ambiguous.csv"),stringsAsFactors=FALSE,check.names=FALSE)
unmapped <- utils::read.csv(file.path(run,"mapping/unmapped.csv"),stringsAsFactors=FALSE,check.names=FALSE)
groups <- complex_groups(expanded,complexes,target)
coverage <- complex_coverage(complexes,groups,target,nonprotein,nested,direct)
frequency <- complex_frequency(coverage,target)
enrichment <- complex_enrichment(coverage,target,background,as.integer(metadata$minimum_overlap),as.numeric(metadata$fdr_cutoff))
navigation <- complex_navigation(groups,coverage,mapping,target)
classes <- c("complete_protein_component_coverage","partial_protein_component_coverage","no_detected_protein_components","not_applicable_no_protein_components")
class_counts <- data.frame(coverage_class=classes,count=as.integer(table(factor(coverage$coverage_class,levels=classes))))
alternatives <- groups[groups$option_count>1,c("complex_id","component_group_id","possible_uniprot_accessions","detected_uniprot_accessions","group_covered"),drop=FALSE]
names(alternatives)[3:4] <- c("possible_accessions","detected_accessions")
stoich <- coverage[,c("complex_id","known_stoichiometry_components","unknown_stoichiometry_components","has_unknown_stoichiometry"),drop=FALSE]
nonprotein_summary <- data.frame(complex_id=as.character(complexes$complex_id),nonprotein_participant_count=integer(nrow(complexes)),participant_types=rep("",nrow(complexes)),stringsAsFactors=FALSE)
nested_summary <- data.frame(complex_id=as.character(complexes$complex_id),nested_complex_count=integer(nrow(complexes)),child_complex_ids=rep("",nrow(complexes)),stringsAsFactors=FALSE)
for(i in seq_len(nrow(complexes))) {
  id <- complexes$complex_id[i];np <- nonprotein[nonprotein$complex_id==id,,drop=FALSE];nc <- nested[nested$parent_complex_id==id,,drop=FALSE]
  nonprotein_summary$nonprotein_participant_count[i] <- nrow(np)
  nonprotein_summary$participant_types[i] <- paste(sort(unique(np$participant_type)),collapse=";")
  nested_summary$nested_complex_count[i] <- nrow(nc)
  nested_summary$child_complex_ids[i] <- paste(sort(unique(nc$child_complex_id)),collapse=";")
}
details <- merge(complexes,coverage,by="complex_id",all.x=TRUE,sort=FALSE)
summary <- data.frame(metric=c("Input rows","Unique canonical proteins","Ambiguous entities","Unmapped entities","Target definition","Target proteins","Background proteins","Complex Portal snapshot","Complex Portal release","Curated complexes","Complexes with protein components","Complexes with detected protein components","Complete protein-component coverage","Partial protein-component coverage","No detected protein components","Not applicable — no protein components","Complexes with alternative groups","Complexes with non-protein participants","Complexes with nested complexes","Complexes with unknown stoichiometry","Complexes represented in target","Complexes tested","Significant complexes","Minimum overlap","FDR method","FDR threshold"),value=as.character(c(nrow(mapping),length(unique(mapping$canonical_uniprot[mapping$mapping_status=="mapped_unique"])),nrow(ambiguous),nrow(unmapped),metadata$target_definition,length(target),length(background),metadata$snapshot_id,metadata$complex_portal_release,nrow(complexes),sum(coverage$total_protein_component_groups>0),sum(coverage$covered_component_groups>0),sum(coverage$coverage_class==classes[1]),sum(coverage$coverage_class==classes[2]),sum(coverage$coverage_class==classes[3]),sum(coverage$coverage_class==classes[4]),sum(coverage$has_alternative_component_groups),sum(coverage$has_nonprotein_participants),sum(coverage$has_nested_complexes),sum(coverage$has_unknown_stoichiometry),sum(frequency$target_member_count>0),nrow(enrichment$all),nrow(enrichment$significant),metadata$minimum_overlap,"Benjamini-Hochberg",metadata$fdr_cutoff)),stringsAsFactors=FALSE)
write_table(groups,"coverage/component_groups.csv")
write_table(alternatives,"coverage/alternative_component_groups.csv")
write_table(coverage,"coverage/complex_coverage.csv")
write_table(class_counts,"coverage/coverage_class_counts.csv")
write_table(direct,"coverage/direct_participants.csv")
write_table(details,"coverage/complex_details.csv")
write_table(stoich,"coverage/stoichiometry_summary.csv")
write_table(nonprotein_summary,"coverage/nonprotein_summary.csv")
write_table(nested_summary,"coverage/nested_complex_summary.csv")
write_table(frequency,"frequency/complex_frequency.csv")
write_table(enrichment$all,"enrichment/complex_enrichment_all.csv")
write_table(enrichment$significant,"enrichment/complex_enrichment_significant.csv")
write_table(enrichment$excluded,"enrichment/complex_enrichment_excluded.csv")
write_table(navigation$protein_to_complexes,"navigation/protein_to_complexes.csv")
write_table(navigation$complex_to_proteins,"navigation/complex_to_proteins.csv")
write_table(summary,"summary.csv")
complex_plots(run,coverage,frequency,enrichment,as.integer(metadata$top_n))
tabs <- list(Summary=summary,`Protein mapping`=mapping,Ambiguous=ambiguous,Unmapped=unmapped,
  `Complex coverage`=coverage,`Coverage classes`=class_counts,`Component groups`=groups,
  `Alternative groups`=alternatives,`Direct participants`=direct,Frequency=frequency,
  Enrichment=enrichment$all,`Enrichment excluded`=enrichment$excluded,
  `Protein to complexes`=navigation$protein_to_complexes,`Complex to proteins`=navigation$complex_to_proteins,
  `Non-protein summary`=nonprotein_summary,`Nested complexes`=nested_summary,Stoichiometry=stoich)
workbook <- openxlsx::createWorkbook()
for(name in names(tabs)){openxlsx::addWorksheet(workbook,name);openxlsx::writeData(workbook,name,tabs[[name]])}
openxlsx::saveWorkbook(workbook,file.path(run,"Complex_analysis.xlsx"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(provenance,"R_session_info.txt"))
cat(jsonlite::toJSON(list(run_id=metadata$run_id,tested=nrow(enrichment$all),significant=nrow(enrichment$significant)),auto_unbox=TRUE),"\n")
