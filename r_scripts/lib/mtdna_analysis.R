# Offline gene-level mtDNA Evidence statistics. No upstream network access.
mtdna_membership <- function(memberships) {
  unique(memberships[,c("entity_key","category"),drop=FALSE])
}

mtdna_frequency <- function(categories,memberships,target,background) {
  membership <- mtdna_membership(memberships)
  rows <- lapply(seq_len(nrow(categories)),function(i) {
    category <- as.character(categories$category[i])
    members <- unique(as.character(membership$entity_key[membership$category==category]))
    hits <- sort(intersect(target,members))
    data.frame(category=category,display_name=as.character(categories$display_name[i]),
      evidence_class=as.character(categories$evidence_class[i]),target_entity_count=length(hits),
      target_size=length(target),fraction_of_target=if(length(target))length(hits)/length(target) else NA_real_,
      background_entity_count=length(intersect(background,members)),background_size=length(background),
      target_entities=paste(hits,collapse=";"),stringsAsFactors=FALSE)
  })
  if(!length(rows))return(data.frame(category=character(),display_name=character(),evidence_class=character(),
    target_entity_count=integer(),target_size=integer(),fraction_of_target=numeric(),
    background_entity_count=integer(),background_size=integer(),target_entities=character()))
  result <- do.call(rbind,rows)
  result[order(-result$target_entity_count,result$category),,drop=FALSE]
}

mtdna_sources <- function(evidence,target) {
  sources <- c("MitoCarta","Gene Ontology","NCBI mtDNA genome")
  frequency <- do.call(rbind,lapply(sources,function(source) {
    hits <- sort(intersect(target,unique(as.character(evidence$entity_key[evidence$source==source]))))
    data.frame(source=source,target_entity_count=length(hits),target_size=length(target),
      fraction_of_target=if(length(target))length(hits)/length(target) else NA_real_,
      entities=paste(hits,collapse=";"),stringsAsFactors=FALSE)
  }))
  agreement <- do.call(rbind,lapply(target,function(key) {
    part <- evidence[evidence$entity_key==key,,drop=FALSE]
    origins <- unique(as.character(part$source));categories <- sort(unique(as.character(part$evidence_category)))
    data.frame(entity_key=key,mitocarta_evidence="MitoCarta" %in% origins,
      go_evidence="Gene Ontology" %in% origins,mtdna_encoded="NCBI mtDNA genome" %in% origins,
      source_count=length(origins),categories=paste(categories,collapse=";"),stringsAsFactors=FALSE)
  }))
  distribution <- data.frame(source_count=0:3,entity_count=as.integer(table(factor(agreement$source_count,levels=0:3))))
  list(frequency=frequency,agreement=agreement,distribution=distribution)
}

mtdna_enrichment <- function(categories,memberships,target,background,minimum_overlap=3L,fdr_cutoff=.05) {
  if(!all(target %in% background))stop("Target outside experimental background")
  membership <- mtdna_membership(memberships)
  reference <- setdiff(background,target)
  columns <- data.frame(category=character(),display_name=character(),Target_count=integer(),
    Target_size=integer(),Reference_count=integer(),Reference_size=integer(),
    Target_nonmember=integer(),Reference_nonmember=integer(),ProteinRatio=numeric(),
    BgRatio=numeric(),p_value=numeric(),FDR=numeric(),target_entities=character())
  excluded_columns <- data.frame(category=character(),display_name=character(),Target_count=integer(),reason=character())
  tested <- list();excluded <- list()
  for(i in seq_len(nrow(categories))) {
    category <- as.character(categories$category[i]);members <- unique(as.character(membership$entity_key[membership$category==category]))
    hits <- sort(intersect(target,members));a <- length(hits)
    if(a<minimum_overlap) {
      excluded[[length(excluded)+1L]] <- data.frame(category=category,
        display_name=as.character(categories$display_name[i]),Target_count=a,reason="below_minimum_overlap")
      next
    }
    c <- length(intersect(reference,members));b <- length(target)-a;d <- length(reference)-c
    p <- stats::fisher.test(matrix(c(a,b,c,d),nrow=2,byrow=TRUE),alternative="greater")$p.value
    tested[[length(tested)+1L]] <- data.frame(category=category,
      display_name=as.character(categories$display_name[i]),Target_count=a,Target_size=length(target),
      Reference_count=c,Reference_size=length(reference),Target_nonmember=b,Reference_nonmember=d,
      ProteinRatio=a/length(target),BgRatio=if(length(reference))c/length(reference) else NA_real_,
      p_value=p,FDR=NA_real_,target_entities=paste(hits,collapse=";"))
  }
  all <- if(length(tested))do.call(rbind,tested) else columns
  if(nrow(all)) {
    all$FDR <- stats::p.adjust(all$p_value,method="BH")
    all <- all[order(all$FDR,all$p_value,all$category),,drop=FALSE]
  }
  list(all=all,significant=all[all$FDR<=fdr_cutoff,,drop=FALSE],
    excluded=if(length(excluded))do.call(rbind,excluded) else excluded_columns)
}

mtdna_comparison <- function(classes,memberships,background) {
  empty <- data.frame(target=character(),category=character(),entity_count=integer(),target_size=integer(),fraction=numeric())
  if(!nrow(classes))return(empty)
  membership <- mtdna_membership(memberships)
  labels <- unique(as.character(classes$classification))
  rows <- list()
  for(label in labels) {
    members <- intersect(background,unique(as.character(classes$entity_key[classes$classification==label])))
    for(category in sort(unique(as.character(membership$category)))) {
      count <- length(intersect(members,membership$entity_key[membership$category==category]))
      rows[[length(rows)+1L]] <- data.frame(target=label,category=category,entity_count=count,
        target_size=length(members),fraction=if(length(members))count/length(members) else NA_real_)
    }
  }
  if(length(rows))do.call(rbind,rows) else empty
}

mtdna_plot <- function(path,draw) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  grDevices::png(paste0(path,".png"),width=1200,height=800,res=140);draw();grDevices::dev.off()
  grDevices::pdf(paste0(path,".pdf"),width=10,height=7);draw();grDevices::dev.off()
}
mtdna_empty <- function(title,message="No observations") {
  graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,message)
}

mtdna_plots <- function(run,frequency,sources,enrichment,evidence,comparison,top_n=20L) {
  at <- function(name)file.path(run,"plots",name)
  mtdna_plot(at("category_frequency"),function() {
    shown <- head(frequency[frequency$target_entity_count>0,,drop=FALSE],top_n)
    if(!nrow(shown))return(mtdna_empty("Category frequency"))
    graphics::barplot(rev(shown$target_entity_count),names.arg=rev(shown$display_name),horiz=TRUE,las=1,
      main="mtDNA evidence category frequency",xlab="Unique target entities",cex.names=.7)
  })
  mtdna_plot(at("source_frequency"),function()graphics::barplot(sources$frequency$target_entity_count,
    names.arg=sources$frequency$source,main="Independent evidence sources",ylab="Unique target entities"))
  mtdna_plot(at("source_count_distribution"),function()graphics::barplot(sources$distribution$entity_count,
    names.arg=as.character(0:3),main="Number of independent evidence sources per target entity",
    xlab="Source count",ylab="Target entities"))
  mtdna_plot(at("enrichment_dot_plot"),function() {
    shown <- enrichment$significant
    title <- "Significant mtDNA evidence category enrichment"
    if(!nrow(shown)) {shown <- enrichment$all;title <- "All tested categories (none met selected FDR)"}
    if(!nrow(shown))return(mtdna_empty(title,"No categories met minimum overlap"))
    graphics::plot(shown$ProteinRatio,seq_len(nrow(shown)),pch=19,
      cex=pmax(1,sqrt(shown$Target_count)),col=grDevices::colorRampPalette(c("red","blue"))(nrow(shown))[rank(shown$FDR)],
      yaxt="n",ylab="Category",xlab="Target entity fraction",main=title)
    graphics::axis(2,at=seq_len(nrow(shown)),labels=shown$display_name,las=1,cex.axis=.7)
  })
  mtdna_plot(at("evidence_records_per_entity"),function() {
    if(!nrow(evidence))return(mtdna_empty("Evidence records per target entity"))
    counts <- sort(table(evidence$entity_key),decreasing=TRUE)
    graphics::barplot(head(counts,top_n),las=2,main="Evidence records per target entity",ylab="Independent records")
  })
  if(length(unique(comparison$target))>=2) mtdna_plot(at("target_category_comparison"),function() {
    matrix <- stats::xtabs(fraction~category+target,data=comparison)
    graphics::image(t(matrix[nrow(matrix):1,,drop=FALSE]),axes=FALSE,main="Descriptive target-category fractions")
    graphics::axis(1,at=seq(0,1,length.out=ncol(matrix)),labels=colnames(matrix),las=2,cex.axis=.7)
    graphics::axis(2,at=seq(0,1,length.out=nrow(matrix)),labels=rev(rownames(matrix)),las=2,cex.axis=.6)
  })
}
