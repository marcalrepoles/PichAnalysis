# Descriptive, offline quantitative proteomics QC. No imputation or exclusion.
qc_transform <- function(raw,transformation) {
  result <- raw
  result[!is.finite(result)] <- NA_real_
  if(transformation=="log2_positive") {
    result[result<=0 & !is.na(result)] <- NA_real_
    result <- log2(result)
  } else if(transformation=="log2p1") {
    if(any(result<0,na.rm=TRUE))stop("Negative values are incompatible with log2p1.")
    result <- log2(result+1)
  } else if(transformation!="none")stop("Unsupported transformation.")
  result
}

qc_detect <- function(raw,zero_is_missing) {
  detected <- is.finite(raw)
  if(zero_is_missing)detected <- detected & raw!=0
  detected[is.na(detected)] <- FALSE
  detected
}

qc_detection <- function(detected,feature,sample) {
  nfeature <- nrow(detected);nsample <- ncol(detected)
  by_sample <- data.frame(sample=sample$sample_id,condition=sample$condition,replicate=sample$replicate,
    total_features=nfeature,detected_features=colSums(detected),
    missing_features=nfeature-colSums(detected),detection_fraction=colMeans(detected),
    missing_fraction=1-colMeans(detected),stringsAsFactors=FALSE)
  by_feature <- data.frame(feature_id=feature$feature_id,display_identifier=feature$display_identifier,
    detected_sample_count=rowSums(detected),sample_count=nsample,
    detection_fraction=rowMeans(detected),stringsAsFactors=FALSE)
  condition_rows <- list()
  for(condition in unique(sample$condition)) {
    columns <- which(sample$condition==condition);count <- rowSums(detected[,columns,drop=FALSE])
    condition_rows[[length(condition_rows)+1L]] <- data.frame(feature_id=feature$feature_id,
      condition=condition,detected_replicates=count,total_replicates=length(columns),
      detection_fraction=count/length(columns),stringsAsFactors=FALSE)
  }
  by_condition <- do.call(rbind,condition_rows)
  consistency <- by_condition
  consistency$classification <- ifelse(consistency$detected_replicates==0,"not_detected",
    ifelse(consistency$detected_replicates==consistency$total_replicates,
      "detected_in_all_replicates","detected_in_some_replicates"))
  list(sample=by_sample,feature=by_feature,condition=by_condition,consistency=consistency)
}

qc_missingness <- function(raw,detected,feature,sample) {
  nfeature <- nrow(raw);nsample <- ncol(raw)
  by_sample <- data.frame(sample=sample$sample_id,missing_count=nfeature-colSums(detected),
    missing_fraction=1-colMeans(detected),zero_count=colSums(raw==0,na.rm=TRUE),
    na_count=colSums(is.na(raw)),nonfinite_count=colSums(is.infinite(raw)),stringsAsFactors=FALSE)
  by_feature <- data.frame(feature_id=feature$feature_id,display_identifier=feature$display_identifier,
    missing_count=nsample-rowSums(detected),missing_fraction=1-rowMeans(detected),
    zero_count=rowSums(raw==0,na.rm=TRUE),na_count=rowSums(is.na(raw)),
    nonfinite_count=rowSums(is.infinite(raw)),stringsAsFactors=FALSE)
  by_condition <- do.call(rbind,lapply(unique(sample$condition),function(condition) {
    part <- by_sample[sample$condition==condition,,drop=FALSE]
    data.frame(condition=condition,sample_count=nrow(part),mean_missing_fraction=mean(part$missing_fraction),
      median_missing_fraction=stats::median(part$missing_fraction),minimum_missing_fraction=min(part$missing_fraction),
      maximum_missing_fraction=max(part$missing_fraction),stringsAsFactors=FALSE)
  }))
  list(sample=by_sample,feature=by_feature,condition=by_condition)
}

qc_distributions <- function(transformed,sample) {
  do.call(rbind,lapply(seq_len(ncol(transformed)),function(i) {
    values <- transformed[,i];values <- values[is.finite(values)]
    quantiles <- if(length(values))as.numeric(stats::quantile(values,c(0,.25,.5,.75,1))) else rep(NA_real_,5)
    data.frame(sample=sample$sample_id[i],condition=sample$condition[i],replicate=sample$replicate[i],
      n_valid=length(values),minimum=quantiles[1],Q1=quantiles[2],median=quantiles[3],
      mean=if(length(values))mean(values) else NA_real_,Q3=quantiles[4],maximum=quantiles[5],
      IQR=if(length(values))stats::IQR(values) else NA_real_,
      SD=if(length(values)>=2)stats::sd(values) else NA_real_,stringsAsFactors=FALSE)
  }))
}

qc_pairwise <- function(transformed,detected,sample) {
  n <- ncol(transformed);pairs <- list();overlaps <- list()
  pearson <- matrix(NA_real_,n,n,dimnames=list(sample$sample_id,sample$sample_id))
  spearman <- pearson;shared <- matrix(0L,n,n,dimnames=dimnames(pearson))
  jaccard <- pearson
  for(i in seq_len(n))for(j in seq_len(n)) {
    valid <- is.finite(transformed[,i]) & is.finite(transformed[,j])
    count <- sum(valid);shared[i,j] <- count
    p <- s <- NA_real_
    if(count>=3 && stats::sd(transformed[valid,i])>0 && stats::sd(transformed[valid,j])>0) {
      p <- suppressWarnings(stats::cor(transformed[valid,i],transformed[valid,j],method="pearson"))
      s <- suppressWarnings(stats::cor(transformed[valid,i],transformed[valid,j],method="spearman"))
      if(!is.finite(s))s <- NA_real_
    }
    pearson[i,j] <- p;spearman[i,j] <- s
    intersection <- sum(detected[,i] & detected[,j]);union <- sum(detected[,i] | detected[,j])
    jaccard[i,j] <- if(union)intersection/union else NA_real_
    if(i<j) {
      pairs[[length(pairs)+1L]] <- data.frame(sample_a=sample$sample_id[i],sample_b=sample$sample_id[j],
        condition_a=sample$condition[i],condition_b=sample$condition[j],n_shared_features=count,
        pearson=p,spearman=s,reason=if(count<3)"fewer_than_three_shared_features" else
          if(is.na(p))"constant_vector" else "",stringsAsFactors=FALSE)
      overlaps[[length(overlaps)+1L]] <- data.frame(sample_a=sample$sample_id[i],sample_b=sample$sample_id[j],
        condition_a=sample$condition[i],condition_b=sample$condition[j],intersection_count=intersection,
        union_count=union,jaccard=jaccard[i,j],stringsAsFactors=FALSE)
    }
  }
  pair_columns <- data.frame(sample_a=character(),sample_b=character(),condition_a=character(),condition_b=character(),
    n_shared_features=integer(),pearson=numeric(),spearman=numeric(),reason=character())
  overlap_columns <- data.frame(sample_a=character(),sample_b=character(),condition_a=character(),condition_b=character(),
    intersection_count=integer(),union_count=integer(),jaccard=numeric())
  pair_frame <- if(length(pairs))do.call(rbind,pairs) else pair_columns
  overlap_frame <- if(length(overlaps))do.call(rbind,overlaps) else overlap_columns
  within <- do.call(rbind,lapply(unique(sample$condition),function(condition) {
    part <- pair_frame[pair_frame$condition_a==condition & pair_frame$condition_b==condition,,drop=FALSE]
    values_p <- part$pearson[is.finite(part$pearson)];values_s <- part$spearman[is.finite(part$spearman)]
    data.frame(condition=condition,replicate_count=sum(sample$condition==condition),pair_count=nrow(part),
      available_pearson_pairs=length(values_p),median_pearson=if(length(values_p))stats::median(values_p) else NA_real_,
      available_spearman_pairs=length(values_s),median_spearman=if(length(values_s))stats::median(values_s) else NA_real_,
      stringsAsFactors=FALSE)
  }))
  list(pairs=pair_frame,pearson=pearson,spearman=spearman,shared=shared,
    overlap=overlap_frame,jaccard=jaccard,within=within)
}

qc_cv <- function(raw,detected,feature,sample) {
  rows <- list()
  for(condition in unique(sample$condition)) {
    indices <- which(sample$condition==condition)
    for(feature_index in seq_len(nrow(raw))) {
      values <- raw[feature_index,indices]
      values <- values[detected[feature_index,indices] & is.finite(values)]
      n <- length(values);average <- if(n)mean(values) else NA_real_
      deviation <- if(n>=2)stats::sd(values) else NA_real_
      cv <- if(n>=2 && average>0)100*deviation/average else NA_real_
      rows[[length(rows)+1L]] <- data.frame(feature_id=feature$feature_id[feature_index],
        display_identifier=feature$display_identifier[feature_index],condition=condition,
        n_observed=n,n_replicates=length(indices),mean=average,sd=deviation,cv_percent=cv,
        stringsAsFactors=FALSE)
    }
  }
  by_feature <- do.call(rbind,rows)
  by_condition <- do.call(rbind,lapply(unique(sample$condition),function(condition) {
    values <- by_feature$cv_percent[by_feature$condition==condition]
    values <- values[is.finite(values)]
    q <- if(length(values))as.numeric(stats::quantile(values,c(.25,.5,.75))) else rep(NA_real_,3)
    data.frame(condition=condition,features_with_cv=length(values),median_cv=q[2],Q1_cv=q[1],
      Q3_cv=q[3],mean_cv=if(length(values))mean(values) else NA_real_,stringsAsFactors=FALSE)
  }))
  list(feature=by_feature,condition=by_condition)
}

qc_pca <- function(transformed,feature,sample) {
  warning <- "Complete-case PCA may emphasize proteins consistently detected across all samples and therefore does not represent proteins with missing quantitative observations."
  complete <- rowSums(is.finite(transformed))==ncol(transformed)
  available <- sum(rowSums(is.finite(transformed))>0)
  matrix <- transformed[complete,,drop=FALSE]
  variable <- if(nrow(matrix))apply(matrix,1,function(row)stats::var(row)>0) else logical()
  variable[is.na(variable)] <- FALSE
  removed <- nrow(matrix)-sum(variable)
  matrix <- matrix[variable,,drop=FALSE]
  status <- if(ncol(transformed)<2)"PCA not available: fewer than two samples" else
    if(nrow(matrix)<2)"PCA not available: fewer than two variable complete-case features" else "Ready"
  empty_scores <- data.frame(sample=sample$sample_id,condition=sample$condition,replicate=sample$replicate)
  empty_loadings <- data.frame(feature_id=character(),PC1=numeric())
  empty_variance <- data.frame(component=character(),variance=numeric(),variance_explained=numeric())
  distances <- matrix(NA_real_,ncol(transformed),ncol(transformed),dimnames=list(sample$sample_id,sample$sample_id))
  if(status=="Ready") {
    analysis <- stats::prcomp(t(matrix),center=TRUE,scale.=FALSE)
    scores <- data.frame(empty_scores,analysis$x,check.names=FALSE)
    loadings <- data.frame(feature_id=feature$feature_id[which(complete)[variable]],analysis$rotation,check.names=FALSE)
    variances <- analysis$sdev^2
    variance <- data.frame(component=names(variances) %||% paste0("PC",seq_along(variances)),
      variance=variances,variance_explained=variances/sum(variances))
    distances <- as.matrix(stats::dist(t(matrix)))
  } else {scores <- empty_scores;loadings <- empty_loadings;variance <- empty_variance}
  summary <- data.frame(status=status,features_available_before_complete_case=available,
    complete_case_features_used=nrow(matrix),zero_variance_features_excluded=removed,
    center=TRUE,scale=FALSE,missingness_policy="complete-case features only; no imputation",
    interpretation_warning=warning,stringsAsFactors=FALSE)
  list(scores=scores,loadings=loadings,variance=variance,summary=summary,distance=distances)
}

`%||%` <- function(left,right)if(is.null(left))right else left

qc_diagnostics <- function(detection,missingness,distribution,pairs,pca,sample) {
  rows <- list();warnings <- list()
  for(i in seq_len(nrow(sample))) {
    id <- sample$sample_id[i];condition <- sample$condition[i]
    within <- pairs$pairs[(pairs$pairs$sample_a==id | pairs$pairs$sample_b==id) &
      pairs$pairs$condition_a==condition & pairs$pairs$condition_b==condition,,drop=FALSE]
    overlap <- pairs$overlap[(pairs$overlap$sample_a==id | pairs$overlap$sample_b==id) &
      pairs$overlap$condition_a==condition & pairs$overlap$condition_b==condition,,drop=FALSE]
    median_or_na <- function(values)if(any(is.finite(values)))stats::median(values[is.finite(values)]) else NA_real_
    score <- pca$scores[pca$scores$sample==id,,drop=FALSE]
    rows[[i]] <- data.frame(sample=id,condition=condition,replicate=sample$replicate[i],
      detected_features=detection$sample$detected_features[i],missing_fraction=missingness$sample$missing_fraction[i],
      median_transformed_value=distribution$median[i],IQR_transformed_value=distribution$IQR[i],
      median_within_condition_pearson=median_or_na(within$pearson),
      median_within_condition_spearman=median_or_na(within$spearman),
      median_within_condition_jaccard=median_or_na(overlap$jaccard),
      PCA_PC1=if("PC1" %in% names(score))score$PC1[1] else NA_real_,
      PCA_PC2=if("PC2" %in% names(score))score$PC2[1] else NA_real_,stringsAsFactors=FALSE)
    add_warning <- function(code,message)warnings[[length(warnings)+1L]] <<- data.frame(
      scope="sample",sample_or_condition=id,warning_code=code,message=message)
    if(detection$sample$detected_features[i]==0)add_warning("no_detected_features","No features were detected under this run's rule.")
    if(!any(is.finite(within$pearson)))add_warning("correlation_unavailable","Within-condition Pearson correlation is unavailable.")
  }
  for(condition in unique(sample$condition))if(sum(sample$condition==condition)<2)
    warnings[[length(warnings)+1L]] <- data.frame(scope="condition",sample_or_condition=condition,
      warning_code="single_replicate_condition",message="Within-condition replicate metrics are unavailable.")
  if(pca$summary$status!="Ready")warnings[[length(warnings)+1L]] <- data.frame(scope="run",sample_or_condition="",
    warning_code="PCA_unavailable",message=pca$summary$status)
  warning_columns <- data.frame(scope=character(),sample_or_condition=character(),warning_code=character(),message=character())
  list(samples=do.call(rbind,rows),warnings=if(length(warnings))do.call(rbind,warnings) else warning_columns)
}

qc_plot <- function(base,draw) {
  dir.create(dirname(base),recursive=TRUE,showWarnings=FALSE)
  grDevices::png(paste0(base,".png"),width=1200,height=850,res=140)
  draw();grDevices::dev.off()
  grDevices::pdf(paste0(base,".pdf"),width=10,height=7)
  draw();grDevices::dev.off()
}
qc_empty_plot <- function(title,message="Not available") {
  graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,message)
}
qc_heatmap <- function(values,title,labels,palette=grDevices::colorRampPalette(c("#2b6cb0","white","#c53030"))(100)) {
  n <- nrow(values)
  graphics::plot.new();graphics::plot.window(xlim=c(0,n),ylim=c(0,n),asp=1)
  finite <- values[is.finite(values)]
  limits <- if(length(finite))range(finite) else c(0,1)
  if(diff(limits)==0)limits <- limits+c(-.5,.5)
  for(i in seq_len(n))for(j in seq_len(n)) {
    value <- values[i,j]
    color <- if(!is.finite(value))"#dddddd" else palette[max(1,min(100,1+floor(99*(value-limits[1])/diff(limits))))]
    graphics::rect(j-1,n-i,j,n-i+1,col=color,border="white")
    graphics::text(j-.5,n-i+.5,if(is.finite(value))sprintf("%.2f",value) else "NA",cex=.65)
  }
  graphics::axis(1,at=seq_len(n)-.5,labels=labels,las=2,cex.axis=.7)
  graphics::axis(2,at=n-(seq_len(n)-.5),labels=labels,las=2,cex.axis=.7)
  graphics::title(main=title)
}
qc_plots <- function(run,raw,transformed,detected,sample,detection,missingness,pairs,cv,pca,transformation) {
  path <- function(name)file.path(run,"plots",name)
  qc_plot(path("detected_features"),function() {
    graphics::barplot(rev(detection$sample$detected_features),names.arg=rev(sample$sample_id),horiz=TRUE,
      xlim=c(0,max(1,nrow(raw))),main="Detected features per sample",xlab="Detected / total original features")
  })
  qc_plot(path("missing_fraction"),function()graphics::barplot(missingness$sample$missing_fraction,
    names.arg=sample$sample_id,ylim=c(0,1),main="Missing fraction per sample",ylab="Fraction"))
  qc_plot(path("sample_distributions"),function() {
    values <- as.vector(transformed);values <- values[is.finite(values)]
    if(!length(values))return(qc_empty_plot(paste("Sample distributions —",transformation)))
    graphics::hist(values,breaks="FD",main=paste("Sample distributions —",transformation),
      xlab="Transformed quantitative value",col="#3182ce")
  })
  qc_plot(path("sample_boxplots"),function() {
    if(!any(is.finite(transformed)))return(qc_empty_plot(paste("Sample boxplots —",transformation)))
    graphics::boxplot(as.data.frame(transformed),names=sample$sample_id,las=2,
      main=paste("Sample boxplots —",transformation),ylab="Transformed quantitative value")
  })
  qc_plot(path("pearson_heatmap"),function()qc_heatmap(pairs$pearson,"Pearson correlation (NA unavailable)",sample$sample_id))
  qc_plot(path("spearman_heatmap"),function()qc_heatmap(pairs$spearman,"Spearman correlation (NA unavailable)",sample$sample_id))
  qc_plot(path("shared_feature_heatmap"),function()qc_heatmap(pairs$shared,"Shared quantitative features",sample$sample_id))
  qc_plot(path("jaccard_heatmap"),function()qc_heatmap(pairs$jaccard,"Detection Jaccard (NA empty union)",sample$sample_id))
  qc_plot(path("pca"),function() {
    if(pca$summary$status!="Ready")return(qc_empty_plot("Complete-case PCA",pca$summary$status))
    scores <- pca$scores;conditions <- unique(scores$condition)
    colors <- seq_along(conditions);names(colors) <- conditions
    x <- scores$PC1;y <- if("PC2" %in% names(scores))scores$PC2 else rep(0,length(x))
    graphics::plot(x,y,col=colors[scores$condition],pch=19,xlab="PC1",ylab=if("PC2" %in% names(scores))"PC2" else "PC2 unavailable",
      main="Complete-case PCA (centered, not unit-variance scaled)")
    graphics::text(x,y,labels=scores$sample,pos=3,cex=.75)
    graphics::legend("topright",legend=conditions,col=colors,pch=19,bty="n")
  })
  qc_plot(path("condition_cv"),function() {
    valid <- cv$feature[is.finite(cv$feature$cv_percent),,drop=FALSE]
    if(!nrow(valid))return(qc_empty_plot("CV by condition","No CV-eligible features"))
    graphics::boxplot(cv_percent~condition,data=valid,main="Feature CV by condition (linear scale)",ylab="CV (%)")
  })
  qc_plot(path("feature_detection_frequency"),function() {
    counts <- table(factor(detection$feature$detected_sample_count,levels=0:ncol(raw)))
    graphics::barplot(counts,names.arg=0:ncol(raw),main="Feature detection frequency",xlab="Samples detected",ylab="Original features")
  })
}