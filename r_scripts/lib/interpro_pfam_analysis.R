split_up <- function(value) {
  if (!length(value) || is.na(value) || !nzchar(trimws(as.character(value)))) return(character())
  values <- unlist(strsplit(trimws(as.character(value)), "[;,|[:space:]]+")); unique(values[nzchar(values) & toupper(values)!="NA"])
}

protein_mapping <- function(catalog, coverage) {
  for(name in c("source_row","original_id")) if(!name %in% names(catalog)) stop("Catalog is missing required columns.")
  if(!"uniprot_accession" %in% names(catalog)) catalog$uniprot_accession<-NA_character_
  rows<-list()
  for(source in unique(catalog$source_row)){
    group<-catalog[catalog$source_row==source,,drop=FALSE];accessions<-unique(unlist(lapply(group$uniprot_accession,split_up)));status<-if(!length(accessions))"unmapped" else if(length(accessions)==1)"mapped_unique" else "ambiguous";candidates<-if(length(accessions))accessions else NA_character_
    for(accession in candidates){annotation<-if(is.na(accession))NA_character_ else coverage$annotation_status[match(accession,coverage$requested_accession)];rows[[length(rows)+1]]<-data.frame(source_row=as.integer(source),original_id=as.character(group$original_id[[1]]),uniprot_accession=accession,canonical_protein_key=if(is.na(accession))NA_character_ else paste0("UP:",accession),mapping_status=status,annotation_status=annotation,source_rows=as.character(source),candidate_accessions=if(status=="ambiguous")paste(accessions,collapse=";") else NA_character_,stringsAsFactors=FALSE)}
  }
  if(!length(rows))return(data.frame(source_row=integer(),original_id=character(),uniprot_accession=character(),canonical_protein_key=character(),mapping_status=character(),annotation_status=character(),source_rows=character(),candidate_accessions=character()))
  mapping<-do.call(rbind,rows);unique_rows<-mapping[mapping$mapping_status=="mapped_unique",,drop=FALSE]
  for(key in unique(unique_rows$canonical_protein_key)){idx<-which(mapping$canonical_protein_key==key & mapping$mapping_status=="mapped_unique");mapping$source_rows[idx]<-paste(sort(unique(mapping$source_row[idx])),collapse=";")}
  mapping
}

unique_proteins <- function(mapping) unique(mapping$uniprot_accession[mapping$mapping_status=="mapped_unique" & !is.na(mapping$uniprot_accession)])
adjust_protein_sets <- function(target,background,authorized=FALSE){target<-unique(target);background<-unique(background);outside<-setdiff(target,background);list(ok=!length(outside)||authorized,target=if(length(outside)&&authorized)intersect(target,background)else target,background=background,outside=outside)}

feature_frequency <- function(membership,target,id,name,type,integrated=NULL){
  pairs<-unique(membership[membership$uniprot_accession%in%target,,drop=FALSE]);columns<-c(id,name,type,if(!is.null(integrated))integrated,"Protein_count","Target_size","Fraction_of_target","Proteins")
  if(!nrow(pairs))return(stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))),columns));groups<-split(pairs$uniprot_accession,pairs[[id]])
  rows<-do.call(rbind,lapply(names(groups),function(feature){part<-pairs[pairs[[id]]==feature,,drop=FALSE];proteins<-sort(unique(part$uniprot_accession));row<-data.frame(ID=feature,Name=part[[name]][[1]],Type=part[[type]][[1]],Protein_count=length(proteins),Target_size=length(target),Fraction_of_target=if(length(target))length(proteins)/length(target) else 0,Proteins=paste(proteins,collapse=";"),stringsAsFactors=FALSE);if(!is.null(integrated))row<-cbind(row[,1:3,drop=FALSE],Integrated=part[[integrated]][[1]],row[,4:ncol(row),drop=FALSE]);row}));names(rows)<-columns;rows[order(-rows$Protein_count,rows[[id]]),,drop=FALSE]
}

feature_enrichment <- function(membership,target,background,id,name,type,minimum=3L,integrated=NULL){
  pairs<-unique(membership[membership$uniprot_accession%in%background,,drop=FALSE]);features<-unique(pairs[[id]]);columns<-c(id,name,type,if(!is.null(integrated))integrated,"Target_count","Target_size","Reference_count","Reference_size","ProteinRatio","BgRatio","p_value","FDR","Proteins")
  empty<-stats::setNames(as.data.frame(matrix(nrow=0,ncol=length(columns))),columns);if(!length(features))return(list(tested=empty,excluded=empty));reference<-setdiff(background,target)
  rows<-do.call(rbind,lapply(features,function(feature){part<-pairs[pairs[[id]]==feature,,drop=FALSE];members<-unique(part$uniprot_accession);hits<-sort(intersect(target,members));refhits<-intersect(reference,members);row<-data.frame(ID=feature,Name=part[[name]][[1]],Type=part[[type]][[1]],Target_count=length(hits),Target_size=length(target),Reference_count=length(refhits),Reference_size=length(reference),ProteinRatio=if(length(target))length(hits)/length(target) else 0,BgRatio=if(length(reference))length(refhits)/length(reference) else 0,p_value=NA_real_,FDR=NA_real_,Proteins=paste(hits,collapse=";"),stringsAsFactors=FALSE);if(!is.null(integrated))row<-cbind(row[,1:3,drop=FALSE],Integrated=part[[integrated]][[1]],row[,4:ncol(row),drop=FALSE]);row}));names(rows)<-columns;eligible<-rows$Target_count>=as.integer(minimum);tested<-rows[eligible,,drop=FALSE];excluded<-rows[!eligible,,drop=FALSE]
  if(nrow(tested)){tested$p_value<-vapply(seq_len(nrow(tested)),function(i){a<-tested$Target_count[i];b<-tested$Target_size[i]-a;c<-tested$Reference_count[i];d<-tested$Reference_size[i]-c;tryCatch(stats::fisher.test(matrix(c(a,c,b,d),nrow=2),alternative="greater")$p.value,error=function(e)NA_real_)},numeric(1));tested$FDR<-stats::p.adjust(tested$p_value,"BH");tested<-tested[order(tested$FDR,tested$p_value,-tested$Target_count),,drop=FALSE]};list(tested=tested,excluded=excluded)
}

location_summary <- function(proteins,interpro,pfam)data.frame(UniProt=proteins,number_of_InterPro_locations=vapply(proteins,function(x)sum(interpro$uniprot_accession==x),integer(1)),number_of_Pfam_locations=vapply(proteins,function(x)sum(pfam$uniprot_accession==x),integer(1)),stringsAsFactors=FALSE)
repeated_features <- function(locations,source,id){if(!nrow(locations))return(data.frame(UniProt=character(),Source_database=character(),Feature_ID=character(),Occurrence_count=integer(),Locations=character()));key<-paste(locations$uniprot_accession,locations[[id]],sep="\r");groups<-split(seq_len(nrow(locations)),key);rows<-lapply(groups,function(ix){if(length(ix)<=1)return(NULL);data.frame(UniProt=locations$uniprot_accession[ix[1]],Source_database=source,Feature_ID=locations[[id]][ix[1]],Occurrence_count=length(ix),Locations=paste(paste0(locations$start[ix],"-",locations$end[ix]),collapse=";"),stringsAsFactors=FALSE)});rows<-Filter(Negate(is.null),rows);if(length(rows))do.call(rbind,rows) else data.frame(UniProt=character(),Source_database=character(),Feature_ID=character(),Occurrence_count=integer(),Locations=character())}

architecture_table <- function(proteins,locations,id,allowed=NULL){
  if(!is.null(allowed)&&"entry_type"%in%names(locations))locations<-locations[tolower(locations$entry_type)%in%allowed,,drop=FALSE]
  rows<-lapply(proteins,function(protein){part<-locations[locations$uniprot_accession==protein,,drop=FALSE];if(nrow(part))part<-part[order(part$start,part$end,part[[id]],part$location_index,part$fragment_index),,drop=FALSE];features<-if(nrow(part))as.character(part[[id]]) else character();overlap<-FALSE;if(nrow(part)>1)overlap<-any(part$start[-1]<=cummax(part$end)[-nrow(part)]);data.frame(UniProt=protein,Architecture=if(length(features))paste(features,collapse=" > ") else NA_character_,Feature_count=length(features),Unique_feature_count=length(unique(features)),Has_repeated_feature=any(duplicated(features)),Has_overlap=overlap,stringsAsFactors=FALSE)});do.call(rbind,rows)
}
architecture_frequency <- function(architecture,target){part<-architecture[architecture$UniProt%in%target & !is.na(architecture$Architecture),,drop=FALSE];if(!nrow(part))return(data.frame(Architecture=character(),Protein_count=integer(),Fraction_of_target=numeric(),Proteins=character()));groups<-split(part$UniProt,part$Architecture);rows<-do.call(rbind,lapply(names(groups),function(name)data.frame(Architecture=name,Protein_count=length(unique(groups[[name]])),Fraction_of_target=length(unique(groups[[name]]))/length(target),Proteins=paste(sort(unique(groups[[name]])),collapse=";"),stringsAsFactors=FALSE)));rows[order(-rows$Protein_count,rows$Architecture),]}
