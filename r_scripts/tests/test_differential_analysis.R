args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",grep("^--file=",args,value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script))),"lib","differential_analysis.R"))
stopifnot(requireNamespace("limma",quietly=TRUE),requireNamespace("statmod",quietly=TRUE))
samples <- data.frame(sample_id=sprintf("S%03d",1:6),condition=rep(c("A","B"),each=3))
design <- da_design(samples,"A","B")
stopifnot(qr(design)$rank==2L,identical(as.numeric(design[,2]),c(.5,.5,.5,-.5,-.5,-.5)))
stopifnot(!grepl("contrasts.fit",paste(deparse(body(da_fit)),collapse=" "),fixed=TRUE))
matrix <- t(vapply(1:25,function(i) {
  baseline <- 10+i
  c(baseline,baseline+1,baseline+3,baseline+(i%%5)-3,
    baseline+(i%%5)-1,baseline+(i%%5)+2)
},numeric(6)))
rownames(matrix) <- paste0("F",1:25);colnames(matrix) <- samples$sample_id
matrix[1,] <- c(10,12,14,7,8,9)
matrix[2,] <- c(10,12,NA,7,8,9)
matrix[3,] <- c(5,5,5,5,5,5)
matrix[4,] <- c(10,12,14,9,11,13)
matrix[5,] <- c(9,11,13,10,12,14)
features <- data.frame(feature_id=rownames(matrix),source_row=1:nrow(matrix),
  display_identifier=paste0("P",1:nrow(matrix)),original_identifier=paste0("P",1:nrow(matrix)))
imputation <- data.frame(feature_id=rownames(matrix),imputed_cells=rep(0L,nrow(matrix)))
fit <- da_fit(matrix,design,FALSE,FALSE)
results <- da_results(matrix,features,samples,fit,"A","B","log2",.05,0,imputation)
stopifnot(abs(results$effect[1]-4)<1e-10,abs(results$effect[2]-3)<1e-10)
stopifnot(results$df_residual[1]==4,results$df_residual[2]==3)
stopifnot(abs(results$log2FC[4]-1)<1e-10,abs(results$fold_change[4]-2)<1e-10,
  abs(results$percent_change[4]-100)<1e-9)
stopifnot(abs(results$log2FC[5]+1)<1e-10,abs(results$fold_change[5]-.5)<1e-10,
  abs(results$percent_change[5]+50)<1e-9)
stopifnot(results$effect_direction[1]=="higher_in_condition_A",
  results$effect_direction[5]=="higher_in_condition_B")
stopifnot(isTRUE(all.equal(results$adj.P.Val[results$tested],
  stats::p.adjust(results$P.Value[results$tested],method="BH"),tolerance=1e-14)))
critical <- stats::qt(.975,fit$df.total[1])
expected_low <- fit$coefficients[1,"A_minus_B"]-
  critical*fit$stdev.unscaled[1,"A_minus_B"]*sqrt(fit$s2.post[1])
stopifnot(abs(results$CI_95_low[1]-expected_low)<1e-10)
stopifnot(all(results$combined_significant==results$fdr_significant))
strict <- da_results(matrix,features,samples,fit,"A","B","log2",.05,10,imputation)
stopifnot(identical(results$adj.P.Val,strict$adj.P.Val),
  !any(strict$combined_significant & !strict$fdr_significant))
continuous <- da_results(matrix,features,samples,fit,"A","B","continuous",.05,0,imputation)
stopifnot(all(is.na(continuous$log2FC)),all(is.na(continuous$fold_change)),
  all(is.na(continuous$percent_change)),abs(continuous$difference_A_minus_B[1]-4)<1e-10)
ranked <- da_rank(results)
stopifnot(nrow(ranked)==nrow(matrix),all(diff(ranked$adj.P.Val[ranked$tested])>=0))
trend <- da_fit(matrix,design,TRUE,FALSE)
robust <- da_fit(matrix,design,FALSE,TRUE)
stopifnot(length(trend$p.value[,2])==nrow(matrix),length(robust$p.value[,2])==nrow(matrix))
# Untested features are retained, excluded from the BH family and assigned no adjusted p-value.
with_untested <- matrix
with_untested[3,] <- c(5,NA,NA,5,NA,NA)
fit_untested <- da_fit(with_untested,design,FALSE,FALSE)
results_untested <- da_results(with_untested,features,samples,fit_untested,
  "A","B","log2",.05,0,imputation)
stopifnot(!results_untested$tested[3],is.na(results_untested$adj.P.Val[3]),
  results_untested$result_class[3]=="Not tested",
  isTRUE(all.equal(results_untested$adj.P.Val[results_untested$tested],
    stats::p.adjust(results_untested$P.Value[results_untested$tested],method="BH"),tolerance=1e-14)))
all_one <- fit
all_one$p.value[,"A_minus_B"] <- 1
all_one$t[,"A_minus_B"] <- 0
one_results <- da_results(matrix,features,samples,all_one,"A","B","log2",.05,0,imputation)
stopifnot(all(one_results$adj.P.Val[one_results$tested]==1),!any(one_results$fdr_significant))
unestimable <- fit
unestimable$coefficients[,"A_minus_B"] <- NA_real_
unestimable$t[,"A_minus_B"] <- NA_real_
unestimable$p.value[,"A_minus_B"] <- NA_real_
stopifnot(inherits(try(da_results(matrix,features,samples,unestimable,
  "A","B","log2",.05,0,imputation),silent=TRUE),"try-error"))
cat("R Differential Statistics tests passed\n")
