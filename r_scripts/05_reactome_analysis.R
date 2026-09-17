args_all <- commandArgs(trailingOnly=FALSE)
script_path <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=", args_all)]), winslash="/", mustWork=TRUE)
script_root <- dirname(script_path)
source(file.path(script_root, "lib", "common.R"))
source(file.path(script_root, "lib", "reactome_analysis.R"))

required_packages <- c("readr", "jsonlite", "ggplot2", "openxlsx")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))

argument_names <- c("catalog", "target-file", "background-file", "output", "project", "pathways",
  "relations", "uniprot-map", "ncbi-map", "diagrams", "summations", "run-id", "snapshot-id",
  "snapshot-manifest", "snapshot-hashes", "reactome-release", "organism", "tax-id", "target-name",
  "background-name", "minimum-overlap", "fdr-cutoff", "top-n", "allow-target-adjustment", "mapping-policy")
p <- parse_arguments(commandArgs(trailingOnly=TRUE), argument_names)
validate_reactome_organism(p$organism, p$`tax-id`)

root <- normalizePath(p$output, winslash="/", mustWork=FALSE)
run_root <- file.path(root, "runs", p$`run-id`)
provenance <- file.path(p$project, "scripts", "runs", paste0(p$`run-id`, "_reactome"))
if (dir.exists(run_root) || dir.exists(provenance)) stop("Reactome run ID already exists.")
for (folder in c("mapping", "pathways", "frequency", "enrichment", "graphs", file.path("raw", p$`run-id`)))
  dir.create(file.path(root, folder), recursive=TRUE, showWarnings=FALSE)
dir.create(run_root, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(provenance, "lib"), recursive=TRUE, showWarnings=FALSE)
file.copy(script_path, file.path(provenance, "05_reactome_analysis.R"), overwrite=FALSE)
file.copy(file.path(script_root, "lib", "reactome_analysis.R"), file.path(provenance, "lib", "reactome_analysis.R"), overwrite=FALSE)
file.copy(file.path(script_root, "lib", "common.R"), file.path(provenance, "lib", "common.R"), overwrite=FALSE)
file.copy(p$`snapshot-manifest`, file.path(provenance, "snapshot_manifest.json"), overwrite=FALSE)
write_json(p, file.path(provenance, "parameters.json"))
write_json(jsonlite::fromJSON(p$`snapshot-hashes`), file.path(provenance, "database_hashes.json"))
writeLines(p$`snapshot-id`, file.path(provenance, "snapshot_id.txt"))
writeLines(p$`reactome-release`, file.path(provenance, "reactome_release.txt"))
capture.output(sessionInfo(), file=file.path(provenance, "session_info.txt"))

catalog <- readr::read_csv(p$catalog, show_col_types=FALSE)
target_catalog <- readr::read_csv(p$`target-file`, show_col_types=FALSE)
background_catalog <- readr::read_csv(p$`background-file`, show_col_types=FALSE)
pathways <- read_reactome_pathways(p$pathways)
uniprot <- read_reactome_mapping(p$`uniprot-map`, "UniProt")
ncbi <- read_reactome_mapping(p$`ncbi-map`, "NCBI Gene fallback")
mapping <- canonical_mapping(catalog, uniprot, ncbi)
target_mapping <- canonical_mapping(target_catalog, uniprot, ncbi)
background_mapping <- canonical_mapping(background_catalog, uniprot, ncbi)
membership <- build_pathway_membership(mapping, uniprot, ncbi, pathways)
relations <- read_reactome_relations(p$relations, pathways)

target_entities_initial <- unique(target_mapping$reactome_entity_key[target_mapping$mapping_status == "mapped_unique"])
background_entities <- unique(background_mapping$reactome_entity_key[background_mapping$mapping_status == "mapped_unique"])
adjustment <- adjust_target_background(target_entities_initial, background_entities,
  identical(tolower(p$`allow-target-adjustment`), "true"))
adjustment_record <- list(
  entities_outside_background=adjustment$outside,
  initial_target_size=length(unique(target_entities_initial)),
  initial_background_size=length(unique(background_entities)),
  final_target_size=length(unique(adjustment$target)),
  final_background_size=length(unique(adjustment$background)),
  entities_removed_or_added=adjustment$outside,
  user_decision=if (length(adjustment$outside)) if (adjustment$ok) "Continue" else "Required" else "Not required")
if (!adjustment$ok) {
  write_json(adjustment_record, file.path(root, "target_background_error.json"))
  write_json(adjustment_record, file.path(root, "raw", p$`run-id`, "target_background_error.json"))
  stop("Target outside background: ", length(adjustment$outside), " canonical Reactome entity/entities are outside the selected background.")
}
unlink(file.path(root, "target_background_error.json"), force=TRUE)
target_entities <- adjustment$target
background_entities <- adjustment$background

frequency <- reactome_frequency(membership, target_entities)
enrichment_result <- reactome_enrichment(membership, target_entities, background_entities, as.integer(p$`minimum-overlap`))
enrichment <- enrichment_result$tested
excluded <- enrichment_result$excluded
significant <- enrichment[!is.na(enrichment$FDR) & enrichment$FDR <= as.numeric(p$`fdr-cutoff`), , drop=FALSE]

mapping_dir <- file.path(root, "mapping")
readr::write_csv(mapping, file.path(mapping_dir, "reactome_entity_mapping.csv"), na="")
readr::write_csv(mapping[mapping$mapping_status == "unmapped", ], file.path(mapping_dir, "unmapped.csv"), na="")
readr::write_csv(mapping[mapping$mapping_status == "ambiguous", ], file.path(mapping_dir, "ambiguous.csv"), na="")
readr::write_csv(membership, file.path(root, "pathways", "pathway_membership.csv"), na="")
readr::write_csv(relations, file.path(root, "pathways", "pathway_hierarchy.csv"), na="")
readr::write_csv(hierarchy_ancestry(relations), file.path(root, "pathways", "pathway_ancestry.csv"), na="")
readr::write_csv(frequency, file.path(root, "frequency", "reactome_frequency.csv"), na="")
readr::write_csv(enrichment, file.path(root, "enrichment", "reactome_enrichment_all.csv"), na="")
readr::write_csv(significant, file.path(root, "enrichment", "reactome_enrichment_significant.csv"), na="")
readr::write_csv(excluded, file.path(root, "enrichment", "enrichment_excluded.csv"), na="")
write_json(adjustment_record, file.path(root, "enrichment", "target_background_adjustment.json"))
readr::write_csv(data.frame(reactome_entity_key=target_entities), file.path(root, "raw", p$`run-id`, "target_entities.csv"))
readr::write_csv(data.frame(reactome_entity_key=background_entities), file.path(root, "raw", p$`run-id`, "background_entities.csv"))

status_by_row <- mapping[!duplicated(mapping$source_row), c("source_row", "mapping_status")]
summary <- data.frame(metric=c(
  "Input entities", "Reactome-mapped entities", "Unmapped entities", "Ambiguous entities",
  "Target definition", "Target size", "Background definition", "Background size",
  "Pathways represented", "Pathways eligible", "Pathways tested", "Significant pathways",
  "Minimum overlap", "FDR threshold", "FDR method", "Reactome release", "Snapshot ID", "Mapping policy",
  "Frequency denominator", "Initial target size", "Initial background size", "Final target size",
  "Final background size", "Entities outside background", "User decision"), value=as.character(c(
  length(unique(catalog$source_row)), sum(status_by_row$mapping_status == "mapped_unique"), sum(status_by_row$mapping_status == "unmapped"),
  sum(status_by_row$mapping_status == "ambiguous"), p$`target-name`, length(target_entities), p$`background-name`,
  length(background_entities), nrow(frequency), sum(frequency$Protein_count >= as.integer(p$`minimum-overlap`)),
  nrow(enrichment), nrow(significant), p$`minimum-overlap`, p$`fdr-cutoff`, "Benjamini-Hochberg",
  p$`reactome-release`, p$`snapshot-id`, p$`mapping-policy`,
  "Unique canonical entities in the selected target", adjustment_record$initial_target_size,
  adjustment_record$initial_background_size, adjustment_record$final_target_size,
  adjustment_record$final_background_size, length(adjustment$outside), adjustment_record$user_decision)),
  stringsAsFactors=FALSE)
readr::write_csv(summary, file.path(root, "summary.csv"))

save_plot <- function(plot, name) {
  ggplot2::ggsave(file.path(root, "graphs", paste0(name, ".png")), plot, width=10, height=7, dpi=160)
  ggplot2::ggsave(file.path(root, "graphs", paste0(name, ".pdf")), plot, width=10, height=7)
}
empty_plot <- function(title, message) ggplot2::ggplot() + ggplot2::annotate("text", x=0, y=0, label=message) +
  ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) + ggplot2::labs(title=title) + ggplot2::theme_void()
top_n <- as.integer(p$`top-n`)
if (nrow(frequency)) {
  display_frequency <- head(frequency[order(-frequency$Protein_count), ], top_n)
  plot <- ggplot2::ggplot(display_frequency, ggplot2::aes(stats::reorder(Pathway, Protein_count), Protein_count)) +
    ggplot2::geom_col(fill="#2563EB") + ggplot2::coord_flip() +
    ggplot2::labs(title="Top Reactome pathways by protein frequency", x="Pathway", y="Unique canonical entities") + ggplot2::theme_minimal()
} else plot <- empty_plot("Top Reactome pathways by protein frequency", "No represented pathways")
save_plot(plot, "reactome_frequency")

display <- if (nrow(significant)) head(significant, top_n) else head(enrichment, top_n)
plot_title <- if (nrow(significant)) "Reactome pathway enrichment" else "Reactome pathways — exploratory display without significance cutoff"
if (nrow(display)) {
  dot <- ggplot2::ggplot(display, ggplot2::aes(GeneRatio, stats::reorder(Pathway, GeneRatio), size=Target_count, color=FDR)) +
    ggplot2::geom_point() + ggplot2::labs(title=plot_title, x="GeneRatio", y="Pathway", color="FDR (BH)") + ggplot2::theme_minimal()
  display$FDR_score <- -log10(pmax(display$FDR, .Machine$double.xmin))
  bar <- ggplot2::ggplot(display, ggplot2::aes(stats::reorder(Pathway, FDR_score), FDR_score)) +
    ggplot2::geom_col(fill="#7C3AED") + ggplot2::coord_flip() +
    ggplot2::labs(title=plot_title, x="Pathway", y="-log10(FDR)") + ggplot2::theme_minimal()
} else {
  dot <- empty_plot(plot_title, "No pathways passed the minimum overlap")
  bar <- empty_plot(plot_title, "No pathways passed the minimum overlap")
}
save_plot(dot, "reactome_enrichment_dot")
save_plot(bar, "reactome_enrichment_fdr")

workbook <- openxlsx::createWorkbook()
tabs <- list(Summary=summary, `Entity mapping`=mapping, Unmapped=mapping[mapping$mapping_status == "unmapped", ],
  Ambiguous=mapping[mapping$mapping_status == "ambiguous", ], `Pathway membership`=membership,
  Hierarchy=relations, Frequency=frequency, `Enrichment all`=enrichment, `Enrichment significant`=significant)
for (name in names(tabs)) {
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, tabs[[name]])
}
openxlsx::saveWorkbook(workbook, file.path(root, "Reactome_analysis.xlsx"), overwrite=TRUE)

metadata <- list(run_id=p$`run-id`, snapshot_id=p$`snapshot-id`, reactome_release=p$`reactome-release`,
  organism=p$organism, tax_id=p$`tax-id`, mapping_policy=p$`mapping-policy`, target_definition=p$`target-name`,
  background_definition=p$`background-name`, minimum_overlap=as.integer(p$`minimum-overlap`),
  fdr_threshold=as.numeric(p$`fdr-cutoff`), fdr_method="Benjamini-Hochberg",
  frequency_definition="Unique canonical entities in the selected target",
  statistical_unit="Unique canonical Reactome entity", source="Local Reactome snapshot",
  network_access=FALSE, target_background_adjustment=adjustment_record)
write_json(metadata, file.path(run_root, "metadata.json"))
write_json(metadata, file.path(root, "latest_metadata.json"))

for (folder in c("mapping", "pathways", "frequency", "enrichment", "graphs")) {
  destination <- file.path(run_root, folder)
  dir.create(destination, recursive=TRUE, showWarnings=FALSE)
  file.copy(list.files(file.path(root, folder), full.names=TRUE), destination, overwrite=TRUE)
}
dir.create(file.path(run_root, "raw"), recursive=TRUE, showWarnings=FALSE)
file.copy(list.files(file.path(root, "raw", p$`run-id`), full.names=TRUE), file.path(run_root, "raw"), overwrite=TRUE)
file.copy(file.path(root, "summary.csv"), file.path(run_root, "summary.csv"), overwrite=TRUE)
file.copy(file.path(root, "Reactome_analysis.xlsx"), file.path(run_root, "Reactome_analysis.xlsx"), overwrite=TRUE)
write_json(metadata, file.path(provenance, "run_metadata.json"))
cat(jsonlite::toJSON(metadata, auto_unbox=TRUE), "\n")
