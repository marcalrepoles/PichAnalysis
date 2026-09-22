string_require <- function() {
  if (!requireNamespace("igraph", quietly=TRUE)) stop("The igraph package is required for STRING analysis.")
  if (!requireNamespace("openxlsx", quietly=TRUE)) stop("The openxlsx package is required for STRING workbook output.")
}
string_graph <- function(nodes, edges) {
  string_require()
  vertices <- data.frame(name=as.character(nodes$string_protein_id), stringsAsFactors=FALSE)
  pairs <- if(nrow(edges)) data.frame(from=as.character(edges$protein_a), to=as.character(edges$protein_b), stringsAsFactors=FALSE) else data.frame(from=character(),to=character())
  igraph::graph_from_data_frame(pairs, directed=FALSE, vertices=vertices)
}
string_node_metrics <- function(nodes, edges) {
  graph <- string_graph(nodes,edges); n <- igraph::vcount(graph)
  if (!n) return(data.frame(string_protein_id=character(),preferred_name=character(),gene_symbol=character(),hop_level=integer(),network_degree=integer(),combined_score_sum=numeric(),combined_score_mean=numeric(),betweenness=numeric(),closeness=numeric(),local_clustering_coefficient=numeric(),component_id=integer(),component_size=integer(),seed_support_count=integer()))
  components <- igraph::components(graph)
  ids <- as.character(nodes$string_protein_id); order <- match(ids,igraph::V(graph)$name)
  degree <- igraph::degree(graph)[order]
  btw <- igraph::betweenness(graph,directed=FALSE,weights=NA)[order]
  cluster <- igraph::transitivity(graph,type="local",isolates="NaN")[order]
  # Closeness is normalized within each connected component. Isolates remain NA.
  close <- rep(NA_real_,n)
  for(component in seq_len(components$no)) {
    ix <- which(components$membership==component)
    if(length(ix)>1) {
      sub <- igraph::induced_subgraph(graph,ix)
      close[match(igraph::V(sub)$name,igraph::V(graph)$name)] <- igraph::closeness(sub,mode="all",weights=NA,normalized=TRUE)
    }
  }
  score_sum <- setNames(rep(0,n),igraph::V(graph)$name)
  for(i in seq_len(nrow(edges))) {
    a <- as.character(edges$protein_a[i]); b <- as.character(edges$protein_b[i]); score <- as.numeric(edges$combined_score[i])
    score_sum[a] <- score_sum[a]+score; score_sum[b] <- score_sum[b]+score
  }
  data.frame(string_protein_id=ids,preferred_name=as.character(nodes$preferred_name),gene_symbol=as.character(nodes$gene_symbol),hop_level=as.integer(nodes$hop_level),network_degree=as.integer(degree),combined_score_sum=as.numeric(score_sum[ids]),combined_score_mean=ifelse(degree>0,as.numeric(score_sum[ids])/degree,NA_real_),betweenness=as.numeric(btw),closeness=as.numeric(close[order]),local_clustering_coefficient=as.numeric(cluster),component_id=as.integer(components$membership[order]),component_size=as.integer(components$csize[components$membership[order]]),seed_support_count=as.integer(nodes$seed_support_count),stringsAsFactors=FALSE)
}
string_components <- function(metrics) {
  if(!nrow(metrics)) return(data.frame(component_id=integer(),component_size=integer(),node_ids=character()))
  parts <- split(metrics$string_protein_id,metrics$component_id)
  data.frame(component_id=as.integer(names(parts)),component_size=as.integer(vapply(parts,length,integer(1))),node_ids=vapply(parts,function(x)paste(sort(x),collapse=";"),character(1)))
}
string_evidence_summary <- function(edges) {
  channels <- setdiff(names(edges),c("protein_a","protein_b","combined_score"))
  if(!length(channels)) return(data.frame(channel=character(),edges_with_nonzero_score=integer(),mean_nonzero_score=numeric(),median_nonzero_score=numeric(),maximum_score=numeric()))
  do.call(rbind,lapply(channels,function(channel) {
    scores <- as.numeric(edges[[channel]]);nonzero <- scores[!is.na(scores)&scores>0]
    data.frame(channel=channel,edges_with_nonzero_score=length(nonzero),mean_nonzero_score=if(length(nonzero))mean(nonzero) else NA_real_,median_nonzero_score=if(length(nonzero))stats::median(nonzero) else NA_real_,maximum_score=if(length(scores))max(scores,na.rm=TRUE) else NA_real_)
  }))
}
string_network_summary <- function(nodes, edges, internal_nodes, internal_edges, metrics, hubs) {
  components <- string_components(metrics); degrees <- metrics$network_degree;n <- nrow(nodes)
  data.frame(seed_count=sum(nodes$hop_level==0),mapped_seed_count=sum(nodes$hop_level==0),internal_node_count=nrow(internal_nodes),internal_edge_count=nrow(internal_edges),degree1_node_count=sum(nodes$hop_level==1),degree2_node_count=sum(nodes$hop_level==2),expanded_node_count=n,expanded_edge_count=nrow(edges),connected_component_count=nrow(components),largest_component_size=if(nrow(components))max(components$component_size) else 0,isolated_seed_count=sum(metrics$network_degree[metrics$hop_level==0]==0),hub_count=nrow(hubs),network_density=if(n>1)2*nrow(edges)/(n*(n-1)) else 0,mean_network_degree=if(length(degrees))mean(degrees) else NA_real_,median_network_degree=if(length(degrees))stats::median(degrees) else NA_real_)
}
string_seed_metrics <- function(seeds, edges, metrics, common) {
  if(!nrow(seeds)) return(data.frame(UniProt=character(),Gene=character(),STRING_ID=character(),direct_neighbor_count=integer(),common_neighbor_count=integer(),network_degree=integer(),betweenness=numeric(),component_id=integer()))
  data.frame(UniProt=as.character(seeds$uniprot_accession),Gene=as.character(seeds$gene_symbol),STRING_ID=as.character(seeds$string_protein_id),direct_neighbor_count=vapply(seeds$string_protein_id,function(id)sum(edges$protein_a==id|edges$protein_b==id),integer(1)),common_neighbor_count=if(nrow(common))rep(nrow(common),nrow(seeds)) else rep(0L,nrow(seeds)),network_degree=metrics$network_degree[match(seeds$string_protein_id,metrics$string_protein_id)],betweenness=metrics$betweenness[match(seeds$string_protein_id,metrics$string_protein_id)],component_id=metrics$component_id[match(seeds$string_protein_id,metrics$string_protein_id)])
}
string_plot <- function(frame, value, label, title, destination, top_n=20L) {
  grDevices::png(paste0(destination,".png"),width=1200,height=800,res=140)
  draw <- function() {
    if(!nrow(frame)) {graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,"No eligible observations");return(invisible(NULL))}
    ordered <- frame[order(frame[[value]],decreasing=TRUE),,drop=FALSE];ordered <- head(ordered,top_n)
    graphics::barplot(rev(ordered[[value]]),names.arg=rev(as.character(ordered[[label]])),horiz=TRUE,las=1,main=title,xlab=value,cex.names=.75,col="#2563EB")
  }
  draw();grDevices::dev.off();grDevices::pdf(paste0(destination,".pdf"),width=10,height=7);draw();grDevices::dev.off()
}
string_distribution_plot <- function(values,title,destination) {
  draw <- function() {if(!length(values)) {graphics::plot.new();graphics::title(main=title);graphics::text(.5,.5,"No eligible observations")} else graphics::barplot(table(values),main=title,xlab="Count",ylab="Nodes",col="#7C3AED")}
  grDevices::png(paste0(destination,".png"),width=1200,height=800,res=140);draw();grDevices::dev.off();grDevices::pdf(paste0(destination,".pdf"),width=10,height=7);draw();grDevices::dev.off()
}
