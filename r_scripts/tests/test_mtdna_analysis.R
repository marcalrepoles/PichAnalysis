script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script <- if(length(script_arg))sub("^--file=","",script_arg[1]) else "r_scripts/tests/test_mtdna_analysis.R"
source(file.path(dirname(script),"..","lib","mtdna_analysis.R"))
categories <- data.frame(category=c("replication","nucleoid","encoded"),display_name=c("Replication","Nucleoid","Encoded"),
  evidence_class=c("replication","localization","genomic_origin"))
membership <- data.frame(entity_key=c("NCBI:1","NCBI:1","NCBI:2","NCBI:3","NCBI:4","NCBI:5"),
  category=c("replication","replication","replication","replication","nucleoid","encoded"))
target <- paste0("NCBI:",1:4)
background <- paste0("NCBI:",1:10)
frequency <- mtdna_frequency(categories,membership,target,background)
stopifnot(frequency$target_entity_count[frequency$category=="replication"]==3L,
  frequency$background_size[1]==10L)
evidence <- data.frame(entity_key=c("NCBI:1","NCBI:1","NCBI:2","NCBI:3","NCBI:4","NCBI:5"),
  source=c("MitoCarta","Gene Ontology","MitoCarta","Gene Ontology","Gene Ontology","NCBI mtDNA genome"),
  evidence_category=c("replication","replication","replication","replication","nucleoid","encoded"))
sources <- mtdna_sources(evidence,target)
stopifnot(sources$agreement$source_count[1]==2L,sum(sources$distribution$entity_count)==4L,
  sources$frequency$target_entity_count[sources$frequency$source=="MitoCarta"]==2L)
enrichment <- mtdna_enrichment(categories,membership,target,background,minimum_overlap=2L,fdr_cutoff=.05)
expected <- fisher.test(matrix(c(3,1,0,6),nrow=2,byrow=TRUE),alternative="greater")$p.value
stopifnot(isTRUE(all.equal(enrichment$all$p_value[1],expected)),
  isTRUE(all.equal(enrichment$all$FDR,p.adjust(enrichment$all$p_value,method="BH"))),
  setequal(enrichment$excluded$category,c("nucleoid","encoded")))
zero <- mtdna_enrichment(categories,membership,target,background,minimum_overlap=10L)
stopifnot(nrow(zero$all)==0L,nrow(zero$significant)==0L,nrow(zero$excluded)==3L)
classes <- data.frame(entity_key=c("NCBI:1","NCBI:2","NCBI:3","NCBI:4"),
  classification=c("A-specific","A-specific","B-specific","B-specific"))
comparison <- mtdna_comparison(classes,membership,background)
stopifnot(comparison$entity_count[comparison$target=="A-specific" & comparison$category=="replication"]==2L,
  comparison$entity_count[comparison$target=="B-specific" & comparison$category=="replication"]==1L,
  !"p_value" %in% names(comparison))
large_background <- paste0("NCBI:",seq_len(100))
large_target <- paste0("NCBI:",seq_len(20))
large_membership <- data.frame(entity_key=paste0("NCBI:",c(1:5,21:25)),category="replication")
large_test <- mtdna_enrichment(categories[1,,drop=FALSE],large_membership,large_target,large_background,3L)
stopifnot(large_test$all$Target_count==5L,large_test$all$Target_size==20L,
  large_test$all$Reference_count==5L,large_test$all$Reference_size==80L)
encoded_only <- mtdna_sources(evidence,"NCBI:5")
stopifnot(encoded_only$agreement$mtdna_encoded,
  !encoded_only$agreement$mitocarta_evidence,!encoded_only$agreement$go_evidence,
  encoded_only$agreement$categories=="encoded")
temporary <- tempfile("mtdna-r-");dir.create(temporary)
mtdna_plots(temporary,frequency,sources,enrichment,evidence,comparison,20L)
stopifnot(file.exists(file.path(temporary,"plots/category_frequency.png")),
  file.exists(file.path(temporary,"plots/category_frequency.pdf")),
  file.exists(file.path(temporary,"plots/target_category_comparison.png")))
cat("R mtDNA Evidence tests passed: dedup, background, sources, Fisher, BH, overlap, zero, comparison, plots\n")
