args <- commandArgs(trailingOnly=FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grepl("^--file=", args)])), ".."), winslash="/")
source(file.path(root, "lib", "reactome_analysis.R"))

tested <- 0L
check <- function(value, label) {
  if (!isTRUE(value)) stop("Reactome R test failed: ", label)
  tested <<- tested + 1L
}
temp <- tempfile("reactome-r-tests-")
dir.create(temp)
writeLines(c(
  "R-HSA-1\tRoot one\tHomo sapiens", "R-HSA-2\tRoot two\tHomo sapiens",
  "R-HSA-3\tShared child\tHomo sapiens", "R-HSA-4\tLeaf\tHomo sapiens",
  "R-MMU-1\tMouse pathway\tMus musculus"), file.path(temp, "ReactomePathways.txt"))
writeLines(c("R-HSA-1\tR-HSA-3", "R-HSA-2\tR-HSA-3", "R-HSA-3\tR-HSA-4", "R-MMU-1\tR-MMU-2"),
  file.path(temp, "ReactomePathwaysRelation.txt"))
writeLines(c(
  "U1\tR-HSA-1\thttps://x\tRoot one\tIEA\tHomo sapiens",
  "U1\tR-HSA-1\thttps://x\tRoot one\tIEA\tHomo sapiens",
  "U1\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens",
  "U2\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens",
  "U3\tR-HSA-1\thttps://x\tRoot one\tIEA\tHomo sapiens",
  "M1\tR-MMU-1\thttps://x\tMouse pathway\tIEA\tMus musculus"), file.path(temp, "UniProt2Reactome.txt"))
writeLines(c(
  "N1\tR-HSA-2\thttps://x\tRoot two\tIEA\tHomo sapiens",
  "N1\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens",
  "N2\tR-HSA-4\thttps://x\tLeaf\tIEA\tHomo sapiens",
  "MN\tR-MMU-1\thttps://x\tMouse pathway\tIEA\tMus musculus"), file.path(temp, "NCBI2Reactome.txt"))

paths <- read_reactome_pathways(file.path(temp, "ReactomePathways.txt"))
up <- read_reactome_mapping(file.path(temp, "UniProt2Reactome.txt"), "UniProt")
ncbi <- read_reactome_mapping(file.path(temp, "NCBI2Reactome.txt"), "NCBI Gene fallback")
check(all(c("U1", "U2", "U3") %in% up$source_id), "1 parser UniProt")
check(all(c("N1", "N2") %in% ncbi$source_id), "2 parser NCBI")
check(all(grepl("^R-HSA-", paths$Reactome_ID)) && !"M1" %in% up$source_id, "3 Homo sapiens filtering")
nonhuman <- tryCatch({validate_reactome_organism("Mus musculus", "10090"); FALSE}, error=function(e) grepl("Unsupported organism", e$message))
check(nonhuman, "4 reject non-human")

catalog <- data.frame(
  source_row=1:7, original_id=c("row1", "row2", "U1;U2", "unknown", "group U1 U2", "row6", "row7"),
  uniprot_accession=c("U1", "", "U1;U2", "", "U1 U2", "U3", ""),
  ncbi_gene_id=c("N1", "N1", "", "BAD", "", "", "N2"),
  gene_symbol=paste0("G", 1:7), protein_name=paste("Protein", 1:7),
  uniprot_function=paste("Function", 1:7), ncbi_summary=paste("Summary", 1:7), stringsAsFactors=FALSE)
mapping <- canonical_mapping(catalog, up, ncbi)
row1 <- mapping[mapping$source_row == 1, ]
check(nrow(row1) == 1 && row1$mapping_route == "UniProt", "5 UniProt precedence")
row2 <- mapping[mapping$source_row == 2, ]
check(row2$mapping_route == "NCBI Gene fallback", "6 NCBI fallback")
check(row1$reactome_entity_key == "UP:U1" && row2$reactome_entity_key == "NCBI:N1", "7 canonical key")
check(nrow(row1) == 1, "8 duplicate route collapse")
check(row1$mapping_status == "mapped_unique", "9 mapped_unique")
check(all(mapping$mapping_status[mapping$source_row == 3] == "ambiguous") && sum(mapping$source_row == 3) == 2, "10 ambiguous")
check(mapping$mapping_status[mapping$source_row == 4] == "unmapped", "11 unmapped")
check(all(mapping$mapping_status[mapping$source_row == 5] == "ambiguous") && sum(mapping$source_row == 5) == 2, "12 protein group ambiguity")
check(row1$uniprot_function == "Function 1" && row1$ncbi_summary == "Summary 1", "12a display annotation preservation")
expanded_group <- rbind(catalog[catalog$source_row == 1, ], catalog[catalog$source_row == 1, ])
expanded_group$uniprot_accession <- c("U1", "U2")
expanded_group$source_row <- 99L
check(all(canonical_mapping(expanded_group, up, ncbi)$mapping_status == "ambiguous"), "12b expanded protein group ambiguity")

membership <- build_pathway_membership(mapping, up, ncbi, paths)
check(!any(duplicated(membership[, c("reactome_entity_key", "Reactome_ID")])), "13 entity/pathway deduplication")
target <- c("UP:U1", "NCBI:N1", "UP:U3")
frequency <- reactome_frequency(membership, target)
check(nrow(frequency) == 3 && frequency$Protein_count[frequency$Reactome_ID == "R-HSA-1"] == 2, "14 frequency")
check(frequency$Protein_fraction[frequency$Reactome_ID == "R-HSA-1"] == 2/3, "15 frequency denominator")

relations <- read_reactome_relations(file.path(temp, "ReactomePathwaysRelation.txt"), paths)
check(nrow(relations) == 3 && all(c("parent_pathway_name", "child_pathway_name") %in% names(relations)), "16 parent child hierarchy")
check(sum(relations$child_pathway_id == "R-HSA-3") == 2, "17 multiple parents")
ancestry <- hierarchy_ancestry(relations)
shared <- ancestry[ancestry$pathway_id == "R-HSA-3", ]
check(shared$top_level_ancestors == "R-HSA-1;R-HSA-2" && shared$minimum_depth == 1, "18 top-level ancestry")
cycle_rel <- data.frame(parent_pathway_id=c("A", "B"), child_pathway_id=c("B", "A"))
check(all(hierarchy_ancestry(cycle_rel)$cycle_detected), "19 cycle protection")

classification <- data.frame(source_row=1:7, classification=c("Shared", "Sporadic", "A-specific", "Not reproducibly detected", "Shared", "Shared", "B-specific"))
check(nrow(select_reactome_target(catalog, "All mapped entities")) == 7, "20 target all mapped")
check(all(select_reactome_target(catalog, "Reproducibly detected", classification)$source_row != 4), "21 presence target")
check(identical(select_reactome_target(catalog, "Manual selection", manual_rows=c(2,7))$source_row, c(2L,7L)), "22 manual target")
default_background <- unique(mapping$reactome_entity_key[mapping$mapping_status == "mapped_unique"])
check(setequal(default_background, c("UP:U1", "NCBI:N1", "UP:U3", "NCBI:N2")), "23 experimental default background")
blocked <- adjust_target_background(c("UP:U1", "UP:OUT"), c("UP:U1"), FALSE)
check(!blocked$ok && identical(blocked$outside, "UP:OUT"), "24 target outside background")
adjusted <- adjust_target_background(c("UP:U1", "UP:OUT"), c("UP:U1"), TRUE)
check(adjusted$ok && identical(adjusted$target, "UP:U1"), "25 authorized adjustment")

enrichment <- reactome_enrichment(membership, target, default_background, minimum_overlap=1)$tested
row <- enrichment[enrichment$Reactome_ID == "R-HSA-1", ]
expected <- phyper(row$Target_count-1, row$Background_count, row$Background_size-row$Background_count, row$Target_size, lower.tail=FALSE)
check(isTRUE(all.equal(row$p_value, expected)), "26 hypergeometric")
minimum <- reactome_enrichment(membership, target, default_background, minimum_overlap=2)
check(all(minimum$tested$Target_count >= 2) && all(minimum$excluded$Target_count < 2), "27 minimum overlap")
check(isTRUE(all.equal(enrichment$FDR, p.adjust(enrichment$p_value, "BH"))), "28 BH correction")
check(nrow(enrichment[enrichment$FDR < 0, , drop=FALSE]) == 0, "29 zero significant valid")
check(all(grepl("^(UP|NCBI):", unlist(strsplit(enrichment$Entities[nzchar(enrichment$Entities)], ";", fixed=TRUE)))), "30 entity list")

cat("R Reactome tests passed:", tested, "checks\n")
