args <- commandArgs(trailingOnly=TRUE)
argument <- function(flag) {i<-match(flag,args);if(is.na(i)||i==length(args))stop(paste("Missing",flag));args[i+1]}
run <- argument("--run"); provenance <- argument("--provenance")
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
script <- if(length(script_arg))sub("^--file=", "", script_arg[1]) else "r_scripts/12_differential_preparation.R"
source(file.path(dirname(script),"lib","differential_preparation.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required for differential preparation.")
if(!requireNamespace("openxlsx",quietly=TRUE))stop("openxlsx is required for differential preparation.")
parameters <- jsonlite::fromJSON(file.path(run,"input","parameters.json"),simplifyVector=TRUE)
samples <- utils::read.csv(file.path(run,"input","sample_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
features <- utils::read.csv(file.path(run,"input","feature_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
input <- utils::read.csv(file.path(run,"input","original_quantitative_matrix.csv"),
  stringsAsFactors=FALSE,check.names=FALSE,colClasses="character")
if(!identical(as.character(input$feature_id),as.character(features$feature_id)))stop("Frozen feature identifiers differ.")
if(!all(samples$sample_id %in% names(input)))stop("Frozen sample metadata differs from matrix.")
raw <- as.matrix(data.frame(lapply(input[samples$sample_id],as.numeric),check.names=FALSE))
rownames(raw) <- features$feature_id; colnames(raw) <- samples$sample_id
observed <- dp_detect(raw,isTRUE(parameters$zero_is_missing))
eligibility <- dp_eligibility(observed,features,samples,parameters$condition_A,parameters$condition_B,
  as.integer(parameters$minimum_observed_per_condition))
if(!any(eligibility$eligibility$continuous_eligible)) {
  dir.create(file.path(run,"eligibility"),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(eligibility$candidates,file.path(run,"eligibility","qualitative_detection_candidates.csv"),row.names=FALSE)
  utils::write.csv(eligibility$eligibility,file.path(run,"eligibility","feature_eligibility.csv"),row.names=FALSE)
  stop("NO_CONTINUOUS_ELIGIBLE: no features are eligible for continuous differential analysis; qualitative candidates are separate.")
}
transformed <- dp_transform(raw,observed,parameters$transformation)
normalization <- dp_normalize(transformed,samples,parameters$normalization)
imputation <- dp_impute(normalization$matrix,eligibility$eligibility$continuous_eligible,samples,
  parameters$imputation_method,parameters$imputation_seed,parameters$imputation_parameters,raw)
eligible_features <- features[eligibility$eligibility$continuous_eligible,,drop=FALSE]
if(anyDuplicated(eligible_features$feature_id) || anyDuplicated(samples$sample_id) ||
   !identical(colnames(imputation$prepared),samples$sample_id) ||
   !identical(rownames(imputation$prepared),eligible_features$feature_id))stop("Prepared matrix identity validation failed.")
if(any(rowSums(is.finite(imputation$prepared))==0) ||
   any(eligibility$eligibility$observed_A[eligibility$eligibility$continuous_eligible]==0) ||
   any(eligibility$eligibility$observed_B[eligibility$eligibility$continuous_eligible]==0))
  stop("Prepared matrix contains an invalid eligible feature.")
contrast <- data.frame(condition_A=parameters$condition_A,condition_B=parameters$condition_B,
  comparison_label=paste(parameters$condition_A,"vs",parameters$condition_B),
  sample_count_A=sum(samples$condition==parameters$condition_A),
  sample_count_B=sum(samples$condition==parameters$condition_B),
  selected_samples=paste(samples$sample_id,collapse=";"),
  comparison_direction="Condition A - Condition B on log scale; Condition A / Condition B on linear scale")
design <- data.frame(sample_id=samples$sample_id,column_name=samples$column_name,
  condition=samples$condition,replicate=samples$replicate,
  Condition_A=as.integer(samples$condition==parameters$condition_A),
  Condition_B=as.integer(samples$condition==parameters$condition_B))
missing_before <- sum(is.na(normalization$matrix[eligibility$eligibility$continuous_eligible,,drop=FALSE]))
metric <- c("Input features","Condition A","Condition B","Samples A","Samples B",
  "Quantification type","Transformation","Normalization","Zero treated as missing",
  "Minimum observed per condition","Input missing cells","Continuous-analysis eligible features",
  "Qualitative detection candidates","Excluded features","Imputation method","Imputation seed",
  "Imputed cells","Remaining missing cells","Prepared matrix features","Prepared matrix samples","Run ID")
value <- c(nrow(raw),parameters$condition_A,parameters$condition_B,contrast$sample_count_A,contrast$sample_count_B,
  parameters$quantification_type,parameters$transformation,parameters$normalization,
  isTRUE(parameters$zero_is_missing),parameters$minimum_observed_per_condition,sum(!observed),
  nrow(eligible_features),nrow(eligibility$candidates),nrow(eligibility$excluded),
  parameters$imputation_method,parameters$imputation_seed,nrow(imputation$cells),
  sum(!is.finite(imputation$prepared)),nrow(imputation$prepared),ncol(imputation$prepared),parameters$run_id)
summary <- data.frame(metric=metric,value=as.character(value),stringsAsFactors=FALSE)
method_parameters <- data.frame(method=parameters$imputation_method,seed=parameters$imputation_seed,
  scope="partial missing values in continuous-analysis eligible features only",
  q=parameters$imputation_parameters$q,sigma=parameters$imputation_parameters$sigma,
  MARGIN=if(parameters$imputation_method=="KNN")1L else parameters$imputation_parameters$MARGIN,
  knn_k_requested=parameters$imputation_parameters$knn_k,
  knn_k_applied=if(parameters$imputation_method=="KNN")
    min(as.integer(parameters$imputation_parameters$knn_k),max(1L,nrow(eligible_features)-1L)) else NA_integer_)
tables <- list(
  "summary.csv"=summary,
  "matrices/original_linear_matrix.csv"=dp_matrix_table(as.matrix(input[samples$sample_id]),features),
  "matrices/transformed_matrix.csv"=dp_matrix_table(transformed,features),
  "matrices/normalized_matrix.csv"=dp_matrix_table(normalization$matrix,features),
  "matrices/prepared_matrix.csv"=dp_matrix_table(imputation$prepared,eligible_features),
  "eligibility/feature_missingness_patterns.csv"=eligibility$patterns,
  "eligibility/feature_eligibility.csv"=eligibility$eligibility,
  "eligibility/qualitative_detection_candidates.csv"=eligibility$candidates,
  "eligibility/excluded_features.csv"=eligibility$excluded,
  "normalization/sample_normalization.csv"=normalization$audit,
  "imputation/imputation_mask.csv"=dp_matrix_table(imputation$mask,eligible_features),
  "imputation/method_parameters.csv"=method_parameters,
  "imputation/imputed_cells.csv"=imputation$cells,
  "imputation/feature_imputation_summary.csv"=imputation$feature_summary,
  "imputation/unresolved_missing_values.csv"=imputation$unresolved,
  "design/contrast_metadata.csv"=contrast,
  "design/design_preview.csv"=design)
for(name in names(tables)) {
  path <- file.path(run,name)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(tables[[name]],path,row.names=FALSE,na="")
}
dp_plots(run,raw,transformed,normalization$matrix,imputation,eligibility,samples,parameters$imputation_method)
workbook <- openxlsx::createWorkbook()
sheets <- list(Summary=summary,`Sample metadata`=samples,`Feature metadata`=features,
  `Feature eligibility`=eligibility$eligibility,`Missingness patterns`=eligibility$patterns,
  `Qualitative candidates`=eligibility$candidates,`Excluded features`=eligibility$excluded,
  Normalization=normalization$audit,`Imputation parameters`=method_parameters,
  `Imputation summary`=imputation$feature_summary,
  `Imputed cells`=imputation$cells,`Unresolved missing`=imputation$unresolved,
  `Contrast metadata`=contrast)
for(name in names(sheets)) {
  openxlsx::addWorksheet(workbook,name)
  openxlsx::writeData(workbook,name,sheets[[name]])
}
openxlsx::saveWorkbook(workbook,file.path(run,"Differential_preparation.xlsx"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(provenance,"R_session_info.txt"))
cat(jsonlite::toJSON(list(run_id=parameters$run_id,eligible=nrow(eligible_features),
  imputed=nrow(imputation$cells)),auto_unbox=TRUE),"\n")
