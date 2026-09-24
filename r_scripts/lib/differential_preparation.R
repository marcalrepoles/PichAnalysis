# Offline matrix preparation only. No differential model or missing-mechanism inference.
dp_detect <- function(raw, zero_is_missing) {
  observed <- is.finite(raw)
  if (zero_is_missing) observed <- observed & raw != 0
  observed[is.na(observed)] <- FALSE
  observed
}

dp_eligibility <- function(observed, features, samples, condition_a, condition_b, minimum) {
  a <- which(samples$condition == condition_a)
  b <- which(samples$condition == condition_b)
  count_a <- rowSums(observed[, a, drop=FALSE])
  count_b <- rowSums(observed[, b, drop=FALSE])
  total_a <- length(a); total_b <- length(b)
  pattern <- character(nrow(observed)); reason <- character(nrow(observed))
  eligible <- logical(nrow(observed)); qualitative <- logical(nrow(observed))
  for (i in seq_len(nrow(observed))) {
    na <- count_a[i]; nb <- count_b[i]
    pattern[i] <- if (na == 0 && nb == 0) "all_missing" else
      if (nb == 0) if (na >= minimum) "detected_only_condition_A" else "sparse_only_condition_A" else
      if (na == 0) if (nb >= minimum) "detected_only_condition_B" else "sparse_only_condition_B" else
      if (na == total_a && nb == total_b) "complete_both_conditions" else
      if (na < minimum && nb < minimum) "insufficient_both_conditions" else "partial_both_conditions"
    eligible[i] <- na >= minimum && nb >= minimum && na + nb > 2
    qualitative[i] <- pattern[i] %in% c("detected_only_condition_A", "detected_only_condition_B")
    reason[i] <- if (eligible[i]) "eligible" else if (na == 0 && nb == 0) "all_missing" else
      if (na == 0) "condition_A_absent" else if (nb == 0) "condition_B_absent" else
      if (na + nb <= 2) "insufficient_total_observations" else
      if (na < minimum && nb < minimum) "insufficient_both_conditions" else
      if (na < minimum) "insufficient_condition_A" else "insufficient_condition_B"
  }
  frame <- data.frame(feature_id=features$feature_id, source_row=features$source_row,
    display_identifier=features$display_identifier, observed_A=count_a, missing_A=total_a-count_a,
    total_A=total_a, observed_B=count_b, missing_B=total_b-count_b, total_B=total_b,
    pattern=pattern, continuous_eligible=eligible, eligibility_reason=reason,
    qualitative_candidate=qualitative, stringsAsFactors=FALSE)
  patterns <- frame[, c("feature_id", "source_row", "display_identifier", "observed_A", "missing_A",
    "total_A", "observed_B", "missing_B", "total_B", "pattern")]
  patterns$detection_fraction_A <- count_a / total_a
  patterns$detection_fraction_B <- count_b / total_b
  candidates <- frame[qualitative, c("feature_id", "display_identifier", "observed_A", "total_A",
    "observed_B", "total_B"), drop=FALSE]
  candidates$condition_detected <- ifelse(frame$observed_A[qualitative] > 0, condition_a, condition_b)
  candidates$classification <- ifelse(frame$observed_A[qualitative] > 0,
    "Detected only in Condition A", "Detected only in Condition B")
  list(eligibility=frame, patterns=patterns, candidates=candidates,
    excluded=frame[!eligible, , drop=FALSE])
}

dp_transform <- function(raw, observed, transformation) {
  transformed <- raw
  transformed[!observed] <- NA_real_
  if (transformation == "log2_positive") {
    if (any(transformed <= 0, na.rm=TRUE)) stop("Nonpositive observed values are incompatible with log2_positive.")
    transformed <- log2(transformed)
  } else if (transformation != "none") stop("Unsupported transformation.")
  transformed
}

dp_normalize <- function(transformed, samples, method) {
  sample_medians <- apply(transformed, 2, function(values) {
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else NA_real_
  })
  reference <- if (any(is.finite(sample_medians))) stats::median(sample_medians[is.finite(sample_medians)]) else NA_real_
  shift <- if (method == "none") rep(0, ncol(transformed)) else
    if (method == "median_center") reference - sample_medians else stop("Unsupported normalization.")
  normalized <- sweep(transformed, 2, ifelse(is.finite(shift), shift, 0), "+")
  medians_after <- apply(normalized, 2, function(values) {
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else NA_real_
  })
  audit <- data.frame(sample=samples$sample_id, condition=samples$condition,
    observed_values=colSums(is.finite(transformed)), median_before=sample_medians,
    reference_median=reference, applied_shift=shift, median_after=medians_after,
    stringsAsFactors=FALSE)
  list(matrix=normalized, audit=audit)
}

dp_impute <- function(normalized, eligible, samples, method, seed, options, raw_original=NULL) {
  input <- normalized[eligible, , drop=FALSE]
  if (!nrow(input)) stop("NO_CONTINUOUS_ELIGIBLE")
  prepared <- input
  if (method != "none") {
    if (!requireNamespace("MsCoreUtils", quietly=TRUE)) stop("IMPUTATION_PACKAGE_UNAVAILABLE: MsCoreUtils")
    if (method %in% c("MinProb", "QRILC") && !requireNamespace("imputeLCMD", quietly=TRUE))
      stop("IMPUTATION_PACKAGE_UNAVAILABLE: imputeLCMD")
    if (method == "KNN" && !requireNamespace("impute", quietly=TRUE))
      stop("IMPUTATION_PACKAGE_UNAVAILABLE: impute")
    if (anyNA(input)) {
      set.seed(as.integer(seed))
    prepared <- tryCatch({
      if (method == "MinProb") MsCoreUtils::impute_matrix(input, method="MinProb",
        q=options$q, sigma=options$sigma, MARGIN=as.integer(options$MARGIN)) else
      if (method == "QRILC") MsCoreUtils::impute_matrix(input, method="QRILC",
        sigma=options$sigma, MARGIN=as.integer(options$MARGIN)) else
      if (method == "KNN") MsCoreUtils::impute_matrix(input, method="knn", MARGIN=1L,
        k=min(as.integer(options$knn_k), max(1L, nrow(input)-1L))) else stop("Unsupported imputation method.")
    }, error=function(error) stop(paste("Imputation failed:", conditionMessage(error))))
    }
    if (!identical(dim(prepared), dim(input))) stop("Imputation changed matrix dimensions.")
    if (!isTRUE(all.equal(as.numeric(prepared[!is.na(input)]), as.numeric(input[!is.na(input)]),
                          tolerance=0))) stop("Imputation modified observed values.")
  }
  rownames(prepared) <- rownames(input); colnames(prepared) <- colnames(input)
  imputed <- is.na(input) & is.finite(prepared)
  unresolved <- is.na(input) & !is.finite(prepared)
  mask <- matrix(0L, nrow(prepared), ncol(prepared), dimnames=dimnames(prepared))
  mask[imputed] <- 1L; mask[unresolved] <- NA_integer_
  cells <- which(imputed, arr.ind=TRUE)
  imputed_cells <- data.frame(feature_id=rownames(prepared)[cells[, 1]],
    sample=colnames(prepared)[cells[, 2]], condition=samples$condition[cells[, 2]],
    original_value=if(is.null(raw_original))rep(NA_real_,nrow(cells)) else raw_original[eligible,,drop=FALSE][cells], prepared_value=prepared[cells],
    imputation_method=rep(method, nrow(cells)), stringsAsFactors=FALSE)
  unresolved_cells <- which(unresolved, arr.ind=TRUE)
  unresolved_table <- data.frame(feature_id=rownames(prepared)[unresolved_cells[, 1]],
    sample=colnames(prepared)[unresolved_cells[, 2]],
    condition=samples$condition[unresolved_cells[, 2]], stringsAsFactors=FALSE)
  feature_summary <- data.frame(feature_id=rownames(prepared), observed_cells=rowSums(is.finite(input)),
    missing_before=rowSums(is.na(input)), imputed_cells=rowSums(imputed),
    missing_after=rowSums(!is.finite(prepared)), stringsAsFactors=FALSE)
  list(prepared=prepared, mask=mask, cells=imputed_cells,
    unresolved=unresolved_table, feature_summary=feature_summary)
}

dp_matrix_table <- function(matrix, features) {
  data.frame(features, as.data.frame(matrix, check.names=FALSE), check.names=FALSE)
}

dp_plot <- function(path, draw) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  grDevices::png(paste0(path,".png"), width=1100, height=750, res=130)
  draw(); grDevices::dev.off()
  grDevices::pdf(paste0(path,".pdf"), width=10, height=7)
  draw(); grDevices::dev.off()
}
dp_empty_plot <- function(title, message) {
  graphics::plot.new(); graphics::title(main=title); graphics::text(.5,.5,message)
}
dp_boxplot <- function(matrix, title) {
  if (!any(is.finite(matrix))) return(dp_empty_plot(title,"No observed values"))
  graphics::boxplot(as.data.frame(matrix), las=2, main=title, ylab="Quantitative value")
}
dp_plots <- function(run, raw, transformed, normalized, imputation, eligibility, samples, method) {
  base <- function(name) file.path(run,"plots",name)
  patterns <- table(eligibility$eligibility$pattern)
  dp_plot(base("missingness_patterns"), function() graphics::barplot(patterns, las=2,
    main="Observed missingness pattern counts", ylab="Original features"))
  dp_plot(base("observed_replicates"), function() {
    graphics::plot(eligibility$eligibility$observed_A, eligibility$eligibility$observed_B,
      xlab="Observed in Condition A", ylab="Observed in Condition B", pch=19,
      main="Observed replicates per condition")
  })
  dp_plot(base("eligibility_counts"), function() graphics::barplot(
    c(Eligible=sum(eligibility$eligibility$continuous_eligible),
      Qualitative=nrow(eligibility$candidates), Excluded=nrow(eligibility$excluded)),
    main="Feature eligibility classification", ylab="Original features"))
  dp_plot(base("before_normalization"), function() dp_boxplot(transformed,"Before normalization"))
  dp_plot(base("after_normalization"), function() dp_boxplot(normalized,"After normalization"))
  dp_plot(base("before_imputation_missingness"), function() graphics::barplot(colSums(is.na(normalized)),
    names.arg=samples$sample_id, main="Missing values before imputation", ylab="Cells"))
  dp_plot(base("imputed_per_sample"), function() graphics::barplot(colSums(imputation$mask==1, na.rm=TRUE),
    names.arg=samples$sample_id, main="Imputed cells per sample", ylab="Cells"))
  dp_plot(base("before_after_imputation"), function() {
    if (method == "none") return(dp_empty_plot("Imputation distribution","No imputation selected"))
    before <- as.numeric(normalized[eligibility$eligibility$continuous_eligible,,drop=FALSE])
    after <- as.numeric(imputation$prepared)
    before <- before[is.finite(before)]; after <- after[is.finite(after)]
    if (!length(before) || !length(after)) return(dp_empty_plot("Imputation distribution","No values available"))
    graphics::boxplot(list(Before=before,After=after), main="Observed distribution before vs after imputation")
  })
  dp_plot(base("qualitative_candidates"), function() {
    counts <- table(factor(eligibility$candidates$classification,
      levels=c("Detected only in Condition A","Detected only in Condition B")))
    graphics::barplot(counts, las=2, main="Qualitative detection candidates", ylab="Original features")
  })
}
