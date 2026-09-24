args <- commandArgs(trailingOnly=TRUE)
argument <- function(flag) {i<-match(flag,args);if(is.na(i)||i==length(args))stop(paste("Missing",flag));args[i+1]}
run <- argument("--run"); provenance <- argument("--provenance")
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script <- if(length(script_arg))sub("^--file=","",script_arg[1]) else "r_scripts/13_differential_analysis.R"
source(file.path(dirname(script),"lib","differential_analysis.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required for differential statistics.")
if(!requireNamespace("openxlsx",quietly=TRUE))stop("openxlsx is required for the differential workbook.")
if(!requireNamespace("limma",quietly=TRUE))stop("LIMMA_UNAVAILABLE")
parameters <- jsonlite::fromJSON(file.path(run,"input","parameters.json"),simplifyVector=TRUE)
samples <- utils::read.csv(file.path(run,"input","sample_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
features <- utils::read.csv(file.path(run,"input","feature_metadata.csv"),stringsAsFactors=FALSE,check.names=FALSE)
input <- utils::read.csv(file.path(run,"input","prepared_matrix.csv"),stringsAsFactors=FALSE,check.names=FALSE)
qualitative <- utils::read.csv(file.path(run,"input","qualitative_detection_candidates.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
imputation_summary <- utils::read.csv(file.path(run,"input","feature_imputation_summary.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
contrast <- utils::read.csv(file.path(run,"input","contrast_metadata.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(input)==0)stop("NO_ESTIMABLE_TESTS: prepared matrix is empty.")
if(anyDuplicated(input$feature_id)||anyDuplicated(samples$sample_id)||
   !identical(tail(names(input),nrow(samples)),samples$sample_id)||
   !all(input$feature_id %in% features$feature_id))stop("Frozen prepared matrix and metadata differ.")
if(nrow(contrast)!=1 || contrast$condition_A!=parameters$condition_A ||
   contrast$condition_B!=parameters$condition_B)stop("Frozen contrast metadata differs.")
if(any(qualitative$feature_id %in% input$feature_id))stop("Qualitative candidates cannot enter the continuous model.")
matrix <- as.matrix(data.frame(lapply(input[samples$sample_id],as.numeric),check.names=FALSE))
rownames(matrix) <- input$feature_id;colnames(matrix) <- samples$sample_id
if(any(is.infinite(matrix)))stop("Prepared matrix contains non-finite values other than NA.")
if(!any(rowSums(is.finite(matrix))>2))stop("NO_ESTIMABLE_TESTS: no feature has residual degrees of freedom.")
design <- da_design(samples,parameters$condition_A,parameters$condition_B)
fitted <- tryCatch(da_fit(matrix,design,isTRUE(parameters$ebayes_trend),
  isTRUE(parameters$ebayes_robust)),error=function(error)stop(paste("Model failure:",conditionMessage(error))))
results <- da_results(matrix,features,samples,fitted,parameters$condition_A,parameters$condition_B,
  parameters$prepared_scale,parameters$fdr_threshold,parameters$minimum_absolute_effect,imputation_summary)
ranked <- da_rank(results)
qualitative$preparation_run_id <- rep(parameters$preparation_run_id,nrow(qualitative))
model_audit <- da_model_audit(results,ncol(matrix))
limma_version <- as.character(utils::packageVersion("limma"))
ebayes_audit <- da_ebayes_audit(results,limma_version,isTRUE(parameters$ebayes_trend),
  isTRUE(parameters$ebayes_robust))
ebayes_audit$r_version <- R.version.string
preparation_metadata <- data.frame(preparation_run_id=parameters$preparation_run_id,
  preparation_run_path=parameters$preparation_run_path,
  preparation_parameters_hash=parameters$preparation_parameters_hash,
  prepared_matrix_hash=parameters$prepared_matrix_hash,
  contrast_metadata_hash=parameters$contrast_metadata_hash,
  transformation=parameters$preparation_transformation,
  normalization=parameters$preparation_normalization,
  imputation=parameters$preparation_imputation,
  imputation_seed=parameters$preparation_imputation_seed,
  minimum_observed_per_condition=parameters$preparation_minimum_observed,
  zero_is_missing=parameters$preparation_zero_is_missing,
  condition_A=parameters$condition_A,condition_B=parameters$condition_B,
  samples_A=parameters$sample_count_A,samples_B=parameters$sample_count_B)
summary <- data.frame(metric=c("Preparation run ID","Condition A","Condition B","Comparison direction",
  "Prepared scale","Samples A","Samples B","Preparation transformation","Preparation normalization",
  "Preparation imputation","Features in prepared matrix","Features containing imputed values",
  "Total imputed cells","Features tested","Features not tested","limma version","eBayes trend",
  "eBayes robust","FDR method","FDR threshold","Minimum absolute effect","Confidence level",
  "FDR-significant features","Combined-significant features","Higher in Condition A",
  "Higher in Condition B","Qualitative detection candidates","Imputation caveat","Run ID"),
  value=as.character(c(parameters$preparation_run_id,parameters$condition_A,parameters$condition_B,
    parameters$comparison_direction,parameters$prepared_scale,parameters$sample_count_A,
    parameters$sample_count_B,parameters$preparation_transformation,parameters$preparation_normalization,
    parameters$preparation_imputation,nrow(results),sum(results$contains_imputed_values),
    sum(results$imputed_cell_count),sum(results$tested),sum(!results$tested),limma_version,
    isTRUE(parameters$ebayes_trend),isTRUE(parameters$ebayes_robust),parameters$fdr_method,
    parameters$fdr_threshold,parameters$minimum_absolute_effect,parameters$confidence_level,
    sum(results$fdr_significant),sum(results$combined_significant),
    sum(results$combined_significant & results$effect>0,na.rm=TRUE),
    sum(results$combined_significant & results$effect<0,na.rm=TRUE),nrow(qualitative),
    parameters$imputation_caveat,parameters$run_id)),stringsAsFactors=FALSE)
tables <- list("summary.csv"=summary,"design/design_matrix.csv"=data.frame(sample=samples$sample_id,design),
  "results/differential_results_all.csv"=results,
  "results/differential_results_tested.csv"=results[results$tested,,drop=FALSE],
  "results/differential_results_untested.csv"=results[!results$tested,,drop=FALSE],
  "results/differential_results_fdr_significant.csv"=results[results$fdr_significant,,drop=FALSE],
  "results/differential_results_significant.csv"=results[results$combined_significant,,drop=FALSE],
  "results/differential_results_ranked.csv"=ranked,
  "qualitative/qualitative_detection_candidates.csv"=qualitative,
  "audit/preparation_run_metadata.csv"=preparation_metadata,
  "audit/model_feature_audit.csv"=model_audit,"audit/ebayes_summary.csv"=ebayes_audit)
for(name in names(tables)) {
  path<-file.path(run,name);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(tables[[name]],path,row.names=FALSE,na="")
}
da_plots(run,results,ranked,parameters)
workbook <- openxlsx::createWorkbook()
sheets <- list(Summary=summary,`Preparation run`=preparation_metadata,
  Design=data.frame(sample=samples$sample_id,design),`All results`=results,
  `Tested results`=results[results$tested,,drop=FALSE],
  `FDR significant`=results[results$fdr_significant,,drop=FALSE],
  Significant=results[results$combined_significant,,drop=FALSE],
  `Not tested`=results[!results$tested,,drop=FALSE],
  `Qualitative candidates`=qualitative,`Model audit`=model_audit,`eBayes audit`=ebayes_audit)
for(name in names(sheets)){openxlsx::addWorksheet(workbook,name);openxlsx::writeData(workbook,name,sheets[[name]])}
openxlsx::saveWorkbook(workbook,file.path(run,"Differential_analysis.xlsx"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(provenance,"R_session_info.txt"))
writeLines(limma_version,file.path(provenance,"limma_version.txt"))
cat(jsonlite::toJSON(list(run_id=parameters$run_id,limma_version=limma_version,
  tested=sum(results$tested),significant=sum(results$combined_significant)),auto_unbox=TRUE),"\n")
