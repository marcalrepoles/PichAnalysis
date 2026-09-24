# Two independent conditions, direct A-minus-B coefficient, no contrasts.fit.
da_design <- function(samples, condition_a, condition_b) {
  if (!identical(unique(samples$condition), c(condition_a, condition_b)))
    stop("Frozen condition order does not match the selected contrast.")
  design <- cbind(Intercept=rep(1, nrow(samples)),
    A_minus_B=ifelse(samples$condition == condition_a, 0.5, -0.5))
  rownames(design) <- samples$sample_id
  if (qr(design)$rank != 2L) stop("DESIGN_NOT_FULL_RANK")
  design
}

da_fit <- function(matrix, design, trend=FALSE, robust=FALSE) {
  if (!requireNamespace("limma",quietly=TRUE)) stop("LIMMA_UNAVAILABLE")
  if (robust && !requireNamespace("statmod",quietly=TRUE)) stop("STATMOD_UNAVAILABLE")
  fit <- limma::lmFit(matrix, design)
  limma::eBayes(fit, trend=trend, robust=robust)
}

da_expand <- function(values, n) {
  if (is.null(values)) return(rep(NA_real_,n))
  if (length(values)==1L) return(rep(as.numeric(values),n))
  as.numeric(values)
}

da_results <- function(matrix, features, samples, fit, condition_a, condition_b,
                       prepared_scale, fdr_threshold, minimum_effect, imputation_summary) {
  n <- nrow(matrix); ids <- rownames(matrix)
  a <- which(samples$condition==condition_a); b <- which(samples$condition==condition_b)
  mean_observed <- function(cols) rowMeans(matrix[,cols,drop=FALSE],na.rm=TRUE)
  mean_a <- mean_observed(a);mean_b <- mean_observed(b)
  effect <- as.numeric(fit$coefficients[,"A_minus_B"])
  t_value <- as.numeric(fit$t[,"A_minus_B"])
  p_value <- as.numeric(fit$p.value[,"A_minus_B"])
  df_residual <- da_expand(fit$df.residual,n)
  df_prior <- da_expand(fit$df.prior,n)
  df_total <- da_expand(fit$df.total,n)
  sigma <- da_expand(fit$sigma,n)
  s2_post <- da_expand(fit$s2.post,n)
  stdev <- as.numeric(fit$stdev.unscaled[,"A_minus_B"])
  b_stat <- if(!is.null(fit$lods))as.numeric(fit$lods[,"A_minus_B"]) else rep(NA_real_,n)
  ave <- da_expand(fit$Amean,n)
  reason <- ifelse(!is.finite(effect),"coefficient_not_estimable",
    ifelse(!is.finite(df_residual) | df_residual<=0,"insufficient_residual_df",
      ifelse(!is.finite(t_value) | !is.finite(p_value),"nonfinite_model_statistic","")))
  tested <- reason==""
  if (!any(tested)) stop("NO_ESTIMABLE_TESTS")
  adjusted <- rep(NA_real_,n)
  adjusted[tested] <- stats::p.adjust(p_value[tested],method="BH")
  standard_error <- stdev*sqrt(s2_post)
  critical <- stats::qt(.975,df=df_total)
  ci_low <- effect-critical*standard_error
  ci_high <- effect+critical*standard_error
  ci_low[!tested] <- NA_real_;ci_high[!tested] <- NA_real_
  imputation_index <- match(ids,imputation_summary$feature_id)
  imputed <- as.integer(imputation_summary$imputed_cells[imputation_index])
  imputed[is.na(imputed)] <- 0L
  safe_fc <- rep(NA_real_,n)
  if (prepared_scale=="log2") {
    possible <- is.finite(effect) & effect<log2(.Machine$double.xmax)
    safe_fc[possible] <- 2^effect[possible]
  }
  fdr_significant <- tested & adjusted <= fdr_threshold
  fdr_significant[is.na(fdr_significant)] <- FALSE
  effect_pass <- tested & abs(effect) >= minimum_effect
  effect_pass[is.na(effect_pass)] <- FALSE
  combined <- fdr_significant & effect_pass
  direction <- ifelse(!is.finite(effect),"not_available",
    ifelse(effect>0,"higher_in_condition_A",ifelse(effect<0,"higher_in_condition_B","no_difference")))
  classification <- ifelse(!tested,"Not tested",
    ifelse(combined & effect>0,paste("Significant — higher in",condition_a),
      ifelse(combined & effect<0,paste("Significant — higher in",condition_b),"Not significant")))
  metadata <- features[match(ids,features$feature_id),,drop=FALSE]
  data.frame(metadata,condition_A=condition_a,condition_B=condition_b,
    n_A_model=rowSums(is.finite(matrix[,a,drop=FALSE])),
    n_B_model=rowSums(is.finite(matrix[,b,drop=FALSE])),
    mean_A_prepared=mean_a,mean_B_prepared=mean_b,
    effect=effect,effect_scale=prepared_scale,
    difference_A_minus_B=if(prepared_scale=="continuous")effect else rep(NA_real_,n),
    log2FC=if(prepared_scale=="log2")effect else rep(NA_real_,n),
    fold_change=safe_fc,percent_change=if(prepared_scale=="log2")(safe_fc-1)*100 else rep(NA_real_,n),
    AveExpr=ave,stdev_unscaled=stdev,sigma=sigma,s2_post=s2_post,
    df_residual=df_residual,df_prior=df_prior,df_total=df_total,
    moderated_t=t_value,P.Value=p_value,adj.P.Val=adjusted,B=b_stat,
    CI_95_low=ci_low,CI_95_high=ci_high,
    imputed_cell_count=imputed,imputed_fraction=imputed/ncol(matrix),
    contains_imputed_values=imputed>0,tested=tested,test_exclusion_reason=reason,
    fdr_significant=fdr_significant,effect_size_pass=effect_pass,
    combined_significant=combined,effect_direction=direction,result_class=classification,
    check.names=FALSE,stringsAsFactors=FALSE)
}

da_rank <- function(results) {
  results[order(!results$tested,ifelse(is.finite(results$adj.P.Val),results$adj.P.Val,Inf),
    ifelse(is.finite(results$P.Value),results$P.Value,Inf),
    ifelse(is.finite(results$effect),-abs(results$effect),Inf),results$feature_id),,drop=FALSE]
}

da_model_audit <- function(results, sample_count) {
  data.frame(feature_id=results$feature_id,n_A_model=results$n_A_model,
    n_B_model=results$n_B_model,
    total_model_observations=results$n_A_model+results$n_B_model,
    missing_model_values=sample_count-results$n_A_model-results$n_B_model,
    imputed_cell_count=results$imputed_cell_count,df_residual=results$df_residual,
    tested=results$tested,test_exclusion_reason=results$test_exclusion_reason)
}

da_ebayes_audit <- function(results, version, trend, robust) {
  finite <- results$df_residual[is.finite(results$df_residual)]
  prior <- results$df_prior[is.finite(results$df_prior)]
  data.frame(limma_version=version,trend=trend,robust=robust,
    features_entered=nrow(results),features_tested=sum(results$tested),
    features_untested=sum(!results$tested),
    residual_df_min=if(length(finite))min(finite) else NA_real_,
    residual_df_median=if(length(finite))stats::median(finite) else NA_real_,
    residual_df_max=if(length(finite))max(finite) else NA_real_,
    df_prior_min=if(length(prior))min(prior) else NA_real_,
    df_prior_median=if(length(prior))stats::median(prior) else NA_real_,
    df_prior_max=if(length(prior))max(prior) else NA_real_)
}

da_plot <- function(path, draw) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  grDevices::png(paste0(path,".png"),width=1100,height=760,res=130)
  draw();grDevices::dev.off()
  grDevices::pdf(paste0(path,".pdf"),width=10,height=7)
  draw();grDevices::dev.off()
}
da_empty <- function(title,message="No values available") {
  graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,message)
}
da_hist <- function(values,title,xlab) {
  values <- values[is.finite(values)]
  if (!length(values)) return(da_empty(title))
  graphics::hist(values,breaks="FD",main=title,xlab=xlab,col="#3182ce")
}

da_plots <- function(run,results,ranked,parameters) {
  path <- function(name)file.path(run,"plots",name)
  tested <- results[results$tested,,drop=FALSE]
  x_label <- if(parameters$prepared_scale=="log2")"log2FC (Condition A - Condition B)" else
    "Difference (Condition A - Condition B)"
  palette <- ifelse(tested$combined_significant & tested$effect>0,"#c53030",
    ifelse(tested$combined_significant & tested$effect<0,"#2b6cb0","#777777"))
  da_plot(path("volcano"),function() {
    if(!nrow(tested))return(da_empty("Effect significance plot"))
    y <- -log10(pmax(tested$adj.P.Val,.Machine$double.xmin))
    graphics::plot(tested$effect,y,col=palette,pch=19,xlab=x_label,ylab="-log10(BH adjusted p-value)",
      main=paste("Volcano —",parameters$condition_A,"vs",parameters$condition_B))
    top <- da_rank(tested)
    top <- head(top,as.integer(parameters$volcano_label_top_n))
    if(nrow(top)) {
      labels <- ifelse(nzchar(top$display_identifier),top$display_identifier,top$feature_id)
      graphics::text(top$effect,-log10(pmax(top$adj.P.Val,.Machine$double.xmin)),labels=labels,pos=3,cex=.65)
    }
  })
  da_plot(path("ma"),function() {
    if(!nrow(tested))return(da_empty("MA plot"))
    graphics::plot(tested$AveExpr,tested$effect,pch=19,col=palette,xlab="Average prepared expression",
      ylab=x_label,main="MA plot");graphics::abline(h=0,lty=2)
  })
  da_plot(path("mean_a_vs_b"),function() {
    if(!nrow(tested))return(da_empty("Mean A vs B"))
    graphics::plot(tested$mean_B_prepared,tested$mean_A_prepared,pch=19,col=palette,
      xlab=paste("Mean",parameters$condition_B),ylab=paste("Mean",parameters$condition_A),
      main="Prepared means: Condition A vs Condition B")
    graphics::abline(a=0,b=1,lty=2)
  })
  da_plot(path("p_value_histogram"),function()da_hist(tested$P.Value,"Raw p-values (tested features)","P.Value"))
  da_plot(path("fdr_histogram"),function()da_hist(tested$adj.P.Val,"BH-adjusted p-values (tested features)","adj.P.Val"))
  da_plot(path("effect_distribution"),function()da_hist(tested$effect,"Differential effect distribution",x_label))
  da_plot(path("residual_df_distribution"),function()da_hist(results$df_residual,"Residual degrees of freedom","df_residual"))
  da_plot(path("top_effects"),function() {
    top <- head(ranked[ranked$tested,,drop=FALSE],as.integer(parameters$top_n))
    if(!nrow(top))return(da_empty("Top differential effects"))
    labels <- ifelse(nzchar(top$display_identifier),top$display_identifier,top$feature_id)
    graphics::barplot(rev(top$effect),names.arg=rev(labels),horiz=TRUE,las=1,
      col=ifelse(rev(top$effect)>0,"#c53030","#2b6cb0"),
      main="Top differential effects",xlab=x_label)
  })
  if(parameters$preparation_imputation!="none")
    da_plot(path("imputation_diagnostic"),function() {
      graphics::plot(results$imputed_fraction,results$effect,pch=19,
        xlab="Fraction of prepared values imputed",ylab=x_label,
        main="Imputed fraction vs differential effect (descriptive)")
    })
}
