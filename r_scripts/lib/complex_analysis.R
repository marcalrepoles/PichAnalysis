# Complex Portal scientific calculations. All inputs are local persisted tables.
complex_options <- function(row) {
  kind <- as.character(row$participant_type)
  raw <- as.character(row$participant_accession_raw)
  canonical <- function(value) {
    # Official Complex Portal exports also name UniProt isoforms and PRO chains.
    # The experimental statistical unit is the canonical UniProt accession.
    if (grepl("^[A-Z0-9]{6}([A-Z0-9]{4})?-(PRO_[0-9]+|[1-9][0-9]*)$", value))
      return(sub("-(PRO_[0-9]+|[1-9][0-9]*)$", "", value))
    value
  }
  if (kind == "protein") return(canonical(as.character(row$uniprot_accession)))
  if (kind == "protein_set" && grepl("^\\[[^]]+\\]$", raw)) {
    options <- trimws(strsplit(sub("^\\[|\\]$", "", raw), ",", fixed=TRUE)[[1]])
    if (!length(options) || any(!nzchar(options))) stop(paste("Invalid Complex Portal component representation:", raw))
    return(sort(unique(vapply(options, canonical, character(1)))))
  }
  if (kind == "other/unknown" && grepl("^[A-Z0-9]{6}([A-Z0-9]{4})?-(PRO_[0-9]+|[1-9][0-9]*)$", raw))
    return(canonical(raw))
  stop(paste("Invalid Complex Portal component representation:", raw))
}

complex_groups <- function(expanded, complexes, target) {
  columns <- data.frame(complex_id=character(),complex_name=character(),component_group_id=character(),
    component_type=character(),source_component_type=character(),source_component_accession_raw=character(),
    component_identity_normalized=logical(),possible_uniprot_accessions=character(),detected_uniprot_accessions=character(),
    group_covered=logical(),option_count=integer(),stoichiometry_raw=character(),
    stoichiometry_known=logical(),occurrence_index=integer(),stringsAsFactors=FALSE)
  if (!nrow(expanded)) return(columns)
  result <- list(); j <- 0L
  for (id in unique(as.character(complexes$complex_id))) {
    part <- expanded[as.character(expanded$complex_id)==id,,drop=FALSE]
    if (!nrow(part)) next
    part <- part[order(as.integer(part$occurrence_index)),,drop=FALSE]
    seen <- character(); group_index <- 0L
    for (i in seq_len(nrow(part))) {
      options <- complex_options(part[i,,drop=FALSE])
      key <- paste(options,collapse=";")
      if (key %in% seen) next
      seen <- c(seen,key);group_index <- group_index+1L
      detected <- sort(intersect(options,target));j <- j+1L
      result[[j]] <- data.frame(complex_id=id,
        complex_name=as.character(complexes$recommended_name[match(id,complexes$complex_id)]),
        component_group_id=sprintf("CG%03d",group_index),
        component_type=if(length(options)>1) "alternative" else "single_protein",
        source_component_type=as.character(part$participant_type[i]),
        source_component_accession_raw=as.character(part$participant_accession_raw[i]),
        component_identity_normalized=as.character(part$participant_type[i]) == "other/unknown" ||
          grepl("-(PRO_[0-9]+|[1-9][0-9]*)$", as.character(part$participant_accession_raw[i])),
        possible_uniprot_accessions=key,detected_uniprot_accessions=paste(detected,collapse=";"),
        group_covered=length(detected)>0L,option_count=length(options),
        stoichiometry_raw=as.character(part$stoichiometry_raw[i]),
        stoichiometry_known=as.character(part$stoichiometry_known[i]) %in% c("1","TRUE","True"),
        occurrence_index=as.integer(part$occurrence_index[i]),stringsAsFactors=FALSE)
    }
  }
  if(!length(result)) return(columns)
  do.call(rbind,result)
}

complex_members <- function(groups, id) {
  subset <- groups[groups$complex_id==id,,drop=FALSE]
  if(!nrow(subset)) return(character())
  sort(unique(unlist(strsplit(paste(subset$possible_uniprot_accessions,collapse=";"),";",fixed=TRUE))))
}

complex_coverage <- function(complexes, groups, target, nonprotein, nested, direct) {
  empty <- data.frame(complex_id=character(),complex_name=character(),total_protein_component_groups=integer(),
    covered_component_groups=integer(),protein_component_coverage=numeric(),coverage_class=character(),
    expected_possible_proteins=character(),detected_member_proteins=character(),
    has_alternative_component_groups=logical(),nonprotein_participant_count=integer(),
    has_nonprotein_participants=logical(),nested_complex_count=integer(),has_nested_complexes=logical(),
    known_stoichiometry_components=integer(),unknown_stoichiometry_components=integer(),
    has_unknown_stoichiometry=logical(),confidence_accession=character(),confidence_name=character(),
    source=character(),stringsAsFactors=FALSE)
  if(!nrow(complexes)) return(empty)
  rows <- lapply(seq_len(nrow(complexes)),function(i) {
    id <- as.character(complexes$complex_id[i]);g <- groups[groups$complex_id==id,,drop=FALSE]
    total <- nrow(g);covered <- sum(g$group_covered)
    class <- if(!total) "not_applicable_no_protein_components" else if(!covered) "no_detected_protein_components" else if(covered==total) "complete_protein_component_coverage" else "partial_protein_component_coverage"
    possible <- complex_members(groups,id);detected <- sort(intersect(possible,target))
    np <- nonprotein[nonprotein$complex_id==id,,drop=FALSE];nest <- nested[nested$parent_complex_id==id,,drop=FALSE]
    dp <- direct[direct$complex_id==id,,drop=FALSE]
    unknown <- sum(as.character(dp$stoichiometry_known) %in% c("0","FALSE","False",""))
    data.frame(complex_id=id,complex_name=as.character(complexes$recommended_name[i]),
      total_protein_component_groups=total,covered_component_groups=covered,
      protein_component_coverage=if(total) covered/total else NA_real_,coverage_class=class,
      expected_possible_proteins=paste(possible,collapse=";"),detected_member_proteins=paste(detected,collapse=";"),
      has_alternative_component_groups=any(g$option_count>1),nonprotein_participant_count=nrow(np),
      has_nonprotein_participants=nrow(np)>0,nested_complex_count=nrow(nest),has_nested_complexes=nrow(nest)>0,
      known_stoichiometry_components=nrow(dp)-unknown,unknown_stoichiometry_components=unknown,
      has_unknown_stoichiometry=unknown>0,confidence_accession=as.character(complexes$confidence_accession[i]),
      confidence_name=as.character(complexes$confidence_name[i]),source=as.character(complexes$source[i]),stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

complex_frequency <- function(coverage,target) {
  result <- data.frame(complex_id=character(),complex_name=character(),target_member_count=integer(),
    target_size=integer(),fraction_of_target=numeric(),protein_component_group_count=integer(),
    possible_member_protein_count=integer(),target_member_proteins=character(),coverage_class=character(),
    protein_component_coverage=numeric(),stringsAsFactors=FALSE)
  if(!nrow(coverage)) return(result)
  do.call(rbind,lapply(seq_len(nrow(coverage)),function(i){r<-coverage[i,,drop=FALSE];members<-if(nzchar(r$expected_possible_proteins))strsplit(r$expected_possible_proteins,";",fixed=TRUE)[[1]] else character();hits<-sort(intersect(members,target));data.frame(complex_id=r$complex_id,complex_name=r$complex_name,target_member_count=length(hits),target_size=length(target),fraction_of_target=if(length(target))length(hits)/length(target) else NA_real_,protein_component_group_count=r$total_protein_component_groups,possible_member_protein_count=length(members),target_member_proteins=paste(hits,collapse=";"),coverage_class=r$coverage_class,protein_component_coverage=r$protein_component_coverage,stringsAsFactors=FALSE)}))
}

complex_enrichment <- function(coverage,target,background,minimum_overlap=3L,fdr_cutoff=.05) {
  if(!all(target %in% background)) stop("Target outside background")
  columns <- data.frame(complex_id=character(),complex_name=character(),Target_count=integer(),Target_size=integer(),
    Reference_count=integer(),Reference_size=integer(),Target_nonmember=integer(),Reference_nonmember=integer(),
    ProteinRatio=numeric(),BgRatio=numeric(),p_value=numeric(),FDR=numeric(),
    target_member_proteins=character(),coverage_class=character(),protein_component_coverage=numeric(),stringsAsFactors=FALSE)
  excluded <- data.frame(complex_id=character(),complex_name=character(),Target_count=integer(),reason=character(),stringsAsFactors=FALSE)
  tested <- list();skipped <- list();a<-0L;b<-0L;reference<-setdiff(background,target)
  for(i in seq_len(nrow(coverage))) {
    row<-coverage[i,,drop=FALSE];members<-if(nzchar(row$expected_possible_proteins))strsplit(row$expected_possible_proteins,";",fixed=TRUE)[[1]] else character()
    hits<-sort(intersect(members,target));target_count<-length(hits)
    reason<-if(!length(members))"no_protein_components" else if(target_count<minimum_overlap)"below_minimum_overlap" else ""
    if(nzchar(reason)){b<-b+1L;skipped[[b]]<-data.frame(complex_id=row$complex_id,complex_name=row$complex_name,Target_count=target_count,reason=reason);next}
    reference_count<-length(intersect(members,reference));target_nonmember<-length(target)-target_count;reference_nonmember<-length(reference)-reference_count
    p<-stats::phyper(target_count-1,target_count+reference_count,target_nonmember+reference_nonmember,length(target),lower.tail=FALSE)
    a<-a+1L;tested[[a]]<-data.frame(complex_id=row$complex_id,complex_name=row$complex_name,Target_count=target_count,Target_size=length(target),Reference_count=reference_count,Reference_size=length(reference),Target_nonmember=target_nonmember,Reference_nonmember=reference_nonmember,ProteinRatio=target_count/length(target),BgRatio=if(length(reference))reference_count/length(reference) else NA_real_,p_value=p,FDR=NA_real_,target_member_proteins=paste(hits,collapse=";"),coverage_class=row$coverage_class,protein_component_coverage=row$protein_component_coverage,stringsAsFactors=FALSE)
  }
  all<-if(length(tested))do.call(rbind,tested) else columns
  if(nrow(all)){all$FDR<-stats::p.adjust(all$p_value,method="BH");all<-all[order(all$FDR,all$p_value,all$complex_id),,drop=FALSE]}
  list(all=all,significant=all[all$FDR<=fdr_cutoff,,drop=FALSE],excluded=if(length(skipped))do.call(rbind,skipped) else excluded)
}

complex_navigation <- function(groups,coverage,mapping,target) {
  p2c <- data.frame(canonical_uniprot=character(),gene_symbol=character(),complex_id=character(),complex_name=character(),component_group_id=character(),alternative_group=logical(),stringsAsFactors=FALSE)
  c2p <- data.frame(complex_id=character(),complex_name=character(),canonical_uniprot=character(),possible_component_protein=logical(),detected_target_protein=logical(),stringsAsFactors=FALSE)
  if(!nrow(groups))return(list(protein_to_complexes=p2c,complex_to_proteins=c2p))
  genes<-setNames(as.character(mapping$gene_symbol),as.character(mapping$canonical_uniprot))
  records<-list();index<-0L
  for(i in seq_len(nrow(groups))) {
    options<-strsplit(groups$possible_uniprot_accessions[i],";",fixed=TRUE)[[1]]
    for(protein in options){index<-index+1L;records[[index]]<-data.frame(canonical_uniprot=protein,gene_symbol=if(!is.na(genes[protein]))genes[protein] else "",complex_id=groups$complex_id[i],complex_name=groups$complex_name[i],component_group_id=groups$component_group_id[i],alternative_group=groups$option_count[i]>1,stringsAsFactors=FALSE)}
  }
  p2c<-unique(do.call(rbind,records));keys<-unique(p2c[,c("complex_id","complex_name","canonical_uniprot")]);c2p<-data.frame(keys,possible_component_protein=TRUE,detected_target_protein=keys$canonical_uniprot %in% target,stringsAsFactors=FALSE)
  list(protein_to_complexes=p2c,complex_to_proteins=c2p)
}

complex_plot <- function(path,draw) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  grDevices::png(paste0(path,".png"),width=1200,height=800,res=140);draw();grDevices::dev.off()
  grDevices::pdf(paste0(path,".pdf"),width=10,height=7);draw();grDevices::dev.off()
}
complex_empty_plot <- function(title,message="No eligible observations") {graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,message)}
complex_plots <- function(run,coverage,frequency,enrichment,top_n=20L) {
  path<-function(name)file.path(run,"plots",name)
  complex_plot(path("coverage_classes"),function(){counts<-table(factor(coverage$coverage_class,levels=c("complete_protein_component_coverage","partial_protein_component_coverage","no_detected_protein_components","not_applicable_no_protein_components")));graphics::barplot(counts,las=2,main="Protein-component coverage classes",ylab="Curated complexes",col="#2563eb",cex.names=.7)})
  complex_plot(path("top_complex_coverage"),function(){x<-coverage[coverage$covered_component_groups>0,,drop=FALSE];if(!nrow(x))return(complex_empty_plot("Top detected complex coverage"));x<-head(x[order(-x$covered_component_groups,-x$protein_component_coverage,-x$total_protein_component_groups,x$complex_id),,drop=FALSE],top_n);labels<-paste0(x$complex_id," (",x$covered_component_groups,"/",x$total_protein_component_groups,")");graphics::barplot(rev(x$covered_component_groups),names.arg=rev(labels),horiz=TRUE,las=1,main="Top detected protein-component groups",xlab="Covered component groups",col="#16a34a",cex.names=.7)})
  complex_plot(path("top_complex_frequency"),function(){x<-frequency[frequency$target_member_count>0,,drop=FALSE];if(!nrow(x))return(complex_empty_plot("Top complex frequency"));x<-head(x[order(-x$target_member_count,x$complex_id),,drop=FALSE],top_n);graphics::barplot(rev(x$target_member_count),names.arg=rev(x$complex_id),horiz=TRUE,las=1,main="Top complex target-member counts",xlab="Unique target proteins",col="#7c3aed",cex.names=.7)})
  complex_plot(path("enrichment_dot_plot"),function(){x<-enrichment$significant;if(!nrow(x))return(complex_empty_plot("Complex enrichment: no significant results"));x<-head(x,top_n);graphics::plot(x$ProteinRatio,seq_len(nrow(x)),pch=21,bg=grDevices::hcl.colors(nrow(x),"YlOrRd",rev=TRUE),cex=1+pmin(x$Target_count,20)/8,yaxt="n",xlab="ProteinRatio",ylab="",main="Significant Complex Portal enrichment (FDR)");graphics::axis(2,at=seq_len(nrow(x)),labels=x$complex_id,las=1)})
  complex_plot(path("complex_size_vs_coverage"),function(){x<-coverage[!is.na(coverage$protein_component_coverage),,drop=FALSE];if(!nrow(x))return(complex_empty_plot("Complex size vs protein-component coverage"));graphics::plot(x$total_protein_component_groups,x$protein_component_coverage,pch=16,col="#2563eb",xlab="Expected protein component groups",ylab="Protein-component coverage",main="Complex size vs coverage")})
  complex_plot(path("target_member_distribution"),function(){if(!nrow(frequency))return(complex_empty_plot("Target-member count distribution"));graphics::barplot(table(frequency$target_member_count),xlab="Unique target member proteins",ylab="Curated complexes",main="Target-member count distribution",col="#059669")})
}
