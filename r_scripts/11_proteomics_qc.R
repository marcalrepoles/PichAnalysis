args <- commandArgs(trailingOnly=TRUE)
argument <- function(flag) {i<-match(flag,args);if(is.na(i)||i==length(args))stop(paste("Missing",flag));args[i+1]}
run <- argument("--run");provenance <- argument("--provenance")
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script <- if(length(script_arg))sub("^--file=","",script_arg[1]) else "r_scripts/11_proteomics_qc.R"
source(file.path(dirname(script),"lib","proteomics_qc.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required for Proteomics QC.")
if(!requireNamespace("openxlsx",quietly=TRUE))stop("openxlsx is required for Proteomics QC workbook output.")
parameters <- jsonlite::fromJSON(file.path(run,"input","parameters.json"),simplifyVector=TRUE)
sample <- utils::read.csv(file.path(run,"input","sample_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
feature <- utils::read.csv(file.path(run,"input","feature_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
input <- utils::read.csv(file.path(run,"input","quantitative_matrix.csv"),stringsAsFactors=FALSE,check.names=FALSE)
if(!identical(as.character(input$feature_id),as.character(feature$feature_id)))stop("Frozen feature identifiers differ.")
if(!all(sample$sample_id %in% names(input)))stop("Frozen sample metadata differs from the matrix.")
raw <- as.matrix(data.frame(lapply(input[sample$sample_id],as.numeric),check.names=FALSE))
colnames(raw) <- sample$sample_id;rownames(raw) <- feature$feature_id
detected <- qc_detect(raw,isTRUE(parameters$zero_is_missing))
transformed <- qc_transform(raw,parameters$transformation)
detection <- qc_detection(detected,feature,sample)
missingness <- qc_missingness(raw,detected,feature,sample)
distribution <- qc_distributions(transformed,sample)
pairs <- qc_pairwise(transformed,detected,sample)
cv <- qc_cv(raw,detected,feature,sample)
pca <- qc_pca(transformed,feature,sample)
diagnostics <- qc_diagnostics(detection,missingness,distribution,pairs,pca,sample)
matrix_table <- function(values)data.frame(feature_id=feature$feature_id,as.data.frame(values,check.names=FALSE),check.names=FALSE)
sample_matrix <- function(values)data.frame(sample=rownames(values),as.data.frame(values,check.names=FALSE),check.names=FALSE)
conditions_with_replicate_correlation <- sum(pairs$within$available_pearson_pairs>0)
summary <- data.frame(metric=c("Input features","Selected samples","Conditions","Quantification type",
  "Transformation","Zero treated as missing","Total quantitative values","Missing quantitative values",
  "Overall missing fraction","Median detected features per sample","Minimum detected features",
  "Maximum detected features","Complete-case features for PCA","PCA available",
  "Conditions with replicate correlation available","Features with at least one CV estimate","Run ID"),
  value=as.character(c(nrow(raw),ncol(raw),length(unique(sample$condition)),parameters$quantification_type,
    parameters$transformation,isTRUE(parameters$zero_is_missing),length(raw),sum(!detected),
    mean(!detected),stats::median(detection$sample$detected_features),
    min(detection$sample$detected_features),max(detection$sample$detected_features),
    pca$summary$complete_case_features_used,pca$summary$status=="Ready",
    conditions_with_replicate_correlation,
    length(unique(cv$feature$feature_id[is.finite(cv$feature$cv_percent)])),parameters$run_id)),
  stringsAsFactors=FALSE)
tables <- list(
  "detection/sample_detection.csv"=detection$sample,
  "detection/feature_detection.csv"=detection$feature,
  "detection/feature_condition_detection.csv"=detection$condition,
  "detection/replicate_consistency.csv"=detection$consistency,
  "missingness/sample_missingness.csv"=missingness$sample,
  "missingness/feature_missingness.csv"=missingness$feature,
  "missingness/condition_missingness.csv"=missingness$condition,
  "distributions/sample_distribution_summary.csv"=distribution,
  "matrices/transformed_matrix.csv"=matrix_table(transformed),
  "matrices/detection_matrix.csv"=matrix_table(1L*detected),
  "replicates/pairwise_correlations.csv"=pairs$pairs,
  "replicates/pearson_correlation_matrix.csv"=sample_matrix(pairs$pearson),
  "replicates/spearman_correlation_matrix.csv"=sample_matrix(pairs$spearman),
  "replicates/shared_feature_count_matrix.csv"=sample_matrix(pairs$shared),
  "replicates/within_condition_correlation_summary.csv"=pairs$within,
  "replicates/detection_overlap.csv"=pairs$overlap,
  "replicates/jaccard_matrix.csv"=sample_matrix(pairs$jaccard),
  "replicates/sample_distance_matrix.csv"=sample_matrix(pca$distance),
  "variability/feature_condition_cv.csv"=cv$feature,
  "variability/condition_cv_summary.csv"=cv$condition,
  "pca/pca_scores.csv"=pca$scores,
  "pca/pca_loadings.csv"=pca$loadings,
  "pca/pca_variance.csv"=pca$variance,
  "pca/pca_summary.csv"=pca$summary,
  "diagnostics/sample_diagnostics.csv"=diagnostics$samples,
  "diagnostics/warnings.csv"=diagnostics$warnings,
  "summary.csv"=summary)
write_table <- function(frame,relative) {
  path<-file.path(run,relative);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(frame,path,row.names=FALSE,na="")
}
for(name in names(tables))write_table(tables[[name]],name)
qc_plots(run,raw,transformed,detected,sample,detection,missingness,pairs,cv,pca,parameters$transformation)
workbook <- openxlsx::createWorkbook()
sheets <- list(Summary=summary,`Sample metadata`=sample,`Sample detection`=detection$sample,
  `Feature detection`=detection$feature,`Condition detection`=detection$condition,
  `Replicate consistency`=detection$consistency,`Sample missingness`=missingness$sample,
  `Feature missingness`=missingness$feature,Distributions=distribution,
  `Pairwise correlations`=pairs$pairs,`Pearson matrix`=sample_matrix(pairs$pearson),
  `Spearman matrix`=sample_matrix(pairs$spearman),`Shared features`=sample_matrix(pairs$shared),
  `Detection overlap`=pairs$overlap,`CV by feature`=cv$feature,`CV summary`=cv$condition,
  `PCA scores`=pca$scores,`PCA variance`=pca$variance,
  `Sample diagnostics`=diagnostics$samples,Warnings=diagnostics$warnings)
for(name in names(sheets)){openxlsx::addWorksheet(workbook,name);openxlsx::writeData(workbook,name,sheets[[name]])}
openxlsx::saveWorkbook(workbook,file.path(run,"Proteomics_QC.xlsx"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(provenance,"R_session_info.txt"))
cat(jsonlite::toJSON(list(run_id=parameters$run_id,features=nrow(raw),samples=ncol(raw),
  pca_status=pca$summary$status),auto_unbox=TRUE),"\n")
