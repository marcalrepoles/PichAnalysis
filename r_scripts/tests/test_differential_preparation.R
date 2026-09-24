args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script))), "lib", "differential_preparation.R"))
samples <- data.frame(sample_id=sprintf("S%03d",1:6),condition=rep(c("A","B"),each=3))
features <- data.frame(feature_id=sprintf("ROW%06d",1:8),source_row=1:8,
  display_identifier=c("P1","P1","P3;P4","P5","P6","P7","P8","P9"))
raw <- rbind(c(10,11,12,20,21,22),c(20,21,22,30,31,32),c(1,2,3,4,5,6),
  c(10,11,NA,8,9,NA),c(10,11,12,NA,NA,NA),c(NA,NA,NA,10,11,12),
  c(10,NA,NA,NA,NA,NA),rep(NA_real_,6))
rownames(raw) <- features$feature_id;colnames(raw) <- samples$sample_id
detected <- dp_detect(raw,TRUE)
classification <- dp_eligibility(detected,features,samples,"A","B",2)
stopifnot(identical(classification$eligibility$pattern,
  c("complete_both_conditions","complete_both_conditions","complete_both_conditions",
    "partial_both_conditions","detected_only_condition_A","detected_only_condition_B",
    "sparse_only_condition_A","all_missing")))
stopifnot(sum(classification$eligibility$continuous_eligible)==4,
  nrow(classification$candidates)==2,
  !any(classification$candidates$feature_id %in% classification$eligibility$feature_id[classification$eligibility$continuous_eligible]))
one <- matrix(c(1,NA,NA,2,NA,NA,1,NA,NA,2,NA,NA),nrow=2)
stopifnot(!any(dp_eligibility(dp_detect(one,TRUE),features[1:2,],samples,"A","B",1)$eligibility$continuous_eligible))
zero <- raw;zero[1,1] <- 0
stopifnot(!dp_detect(zero,TRUE)[1,1],dp_detect(zero,FALSE)[1,1])
stopifnot(is.na(dp_transform(zero,dp_detect(zero,TRUE),"log2_positive")[1,1]))
stopifnot(inherits(try(dp_transform(zero,dp_detect(zero,FALSE),"log2_positive"),silent=TRUE),"try-error"))
transformed <- dp_transform(raw,detected,"log2_positive")
normal_none <- dp_normalize(transformed,samples,"none")
stopifnot(identical(normal_none$matrix,transformed))
simple <- matrix(c(1,3,5,7,NA,9),nrow=3,dimnames=list(NULL,c("S001","S002")))
simple_samples <- samples[1:2,,drop=FALSE]
median_center <- dp_normalize(simple,simple_samples,"median_center")
stopifnot(abs(median_center$audit$applied_shift[1]-2.5)<1e-12,
  abs(median_center$audit$applied_shift[2]+2.5)<1e-12,
  is.na(median_center$matrix[2,2]))
no_imp <- dp_impute(normal_none$matrix,classification$eligibility$continuous_eligible,
  samples,"none",12345,list(q=.01,sigma=1,MARGIN=2,knn_k=3))
stopifnot(anyNA(no_imp$prepared),nrow(no_imp$cells)==0,
  !"ROW000005" %in% rownames(no_imp$prepared))

# Exercise all three official MsCoreUtils methods on a sufficiently broad fixture.
stopifnot(requireNamespace("MsCoreUtils",quietly=TRUE),
  requireNamespace("imputeLCMD",quietly=TRUE),requireNamespace("impute",quietly=TRUE))
wide <- matrix(log2(seq(10,129)),nrow=20,ncol=6,
  dimnames=list(paste0("F",1:20),samples$sample_id))
wide[1,3] <- NA_real_;wide[2,6] <- NA_real_;wide[3,1] <- NA_real_
options <- list(q=.01,sigma=1,MARGIN=2,knn_k=3)
for(method in c("MinProb","QRILC","KNN")) {
  result <- dp_impute(wide,rep(TRUE,nrow(wide)),samples,method,12345,options)
  again <- dp_impute(wide,rep(TRUE,nrow(wide)),samples,method,12345,options)
  stopifnot(all(is.finite(result$prepared)),nrow(result$cells)==3,
    isTRUE(all.equal(result$prepared,again$prepared,tolerance=0)),
    isTRUE(all.equal(as.numeric(result$prepared[!is.na(wide)]),as.numeric(wide[!is.na(wide)]),tolerance=0)))
}
# A 1+2 observed feature meets the residual requirement at explicit minimum=1.
small <- matrix(c(1,NA,NA,2,3,NA),nrow=1,
  dimnames=list("F1",samples$sample_id))
small_feature <- data.frame(feature_id="F1",source_row=1,display_identifier="F1")
stopifnot(dp_eligibility(dp_detect(small,TRUE),small_feature,samples,"A","B",1)$eligibility$continuous_eligible)
first_seed <- dp_impute(wide,rep(TRUE,nrow(wide)),samples,"MinProb",12345,options)
second_seed <- dp_impute(wide,rep(TRUE,nrow(wide)),samples,"MinProb",54321,options)
stopifnot(!isTRUE(all.equal(first_seed$prepared[is.na(wide)],second_seed$prepared[is.na(wide)])),
  isTRUE(all.equal(first_seed$prepared[!is.na(wide)],second_seed$prepared[!is.na(wide)],tolerance=0)))
cat("R Differential Preparation tests passed\n")
