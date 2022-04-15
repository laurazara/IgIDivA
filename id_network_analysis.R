
# load libraries ----------------------------------------------------------------

library(data.table)
library(stringr)
library(dplyr)
library(readxl)
library(ggsci)
library(tidygraph)
library(ggraph)
library(igraph)
source('main_block_selection.R')
source('amino_acids_mutation.R')

# read input files ----------------------------------------------------------------

doGraph <- function(highly_sim_clonos_file,grouped_alignment_file,sample_id,save_path=getwd(),include_jump=TRUE, col_start = 5, col_end = 313, min_reads = 10, highly_sim_clonos = c(1), nodes_size_scaling = TRUE, include_aa_muts = TRUE){
  #for FR1 samples
  #col_start=59
  cys_pre_cdr3_length = 3
  col_pos1 = 5 #start position to begin directly from
  #the sequence alignment
  
  save_path = paste0(save_path,'/',sample_id)
  dir.create(save_path,showWarnings = FALSE)
  
  high_sim = data.table::fread(highly_sim_clonos_file,
                               header = TRUE,
                               sep = "\t",
                               stringsAsFactors = FALSE)
  
  alignment = data.table::fread(grouped_alignment_file, 
                                header = TRUE, 
                                sep = "\t", 
                                stringsAsFactors = FALSE)
  #colnames(alignment)[col_start:col_end] = as.character(1:(col_end-col_start+1))
  #for FR1
  colnames(alignment)[col_start:col_end] = as.character(colnames(alignment)[col_start:col_end])
  # find the germline ---------------------------------------------------------------
  gene_count = alignment[,.(N=sum(N)),by=V.GENE.and.allele]
  gene_max_N = gene_count$V.GENE.and.allele[which(gene_count$N==max(gene_count$N))]
  germline = alignment[which(alignment$V.GENE.and.allele==gene_max_N), ]
  germline = germline[which(germline$cluster_id == "-"), ]
  
  
  # pre processing ------------------------------------------------------------------
  
  related_clonos = high_sim[highly_sim_clonos, ]$prev_cluster
  related_clonos = stringr::str_squish(related_clonos)
  related_clonos = stringr::str_split(related_clonos, " ")
  related_clonos_final = character(0)
  for (i in 1:length(highly_sim_clonos)) {
    related_clonos_temp = related_clonos[[i]]
    related_clonos_temp = c(related_clonos_temp)
    related_clonos_final = unique(append(related_clonos_final, related_clonos_temp))
  }
  alignment = alignment[which(alignment$cluster_id %in% related_clonos_final), ]
  joined_filtered = alignment[which(alignment$V.GENE.and.allele == germline$V.GENE.and.allele), ]
  
  # General calculations
  number_related_clonos = length(related_clonos_final)
  total_nt_variants = nrow(joined_filtered)
  total_seqs = sum(joined_filtered$N)
  singletons =  nrow(joined_filtered[which(joined_filtered$N == 1), ])
  expanded_nt_vars = nrow(joined_filtered[which(joined_filtered$N != 1), ])
  nonsingl = (joined_filtered[which(joined_filtered$N != 1), ])
  expanded_seqs = sum(nonsingl$N)
  N_main_nt_var = max(joined_filtered$N)
  summary = data.table("Number_related_clonos" = number_related_clonos, "Number_nt_vars" = total_nt_variants, "Total_seqs"= total_seqs, "Singletons" = singletons, "Expanded_nt_vars" = expanded_nt_vars, "Expanded_seqs" = expanded_seqs, "N_main_nt_var" = N_main_nt_var)
  write.table(summary, paste0(save_path,"/summary-calculations_",sample_id,".txt"), sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  #Select main nt var and its muts ---------------------------------------------------------------
  main = joined_filtered[which.max(joined_filtered$N), ]
  main = main[, col_start:col_end]
  main = as.data.table(t(main))
  #main$pos = as.character(1:nrow(main))
  #For FR1
  main$pos = as.character(colnames(alignment)[col_start:col_end])
  main_gaps = main[which(main$V1 == "."), ] 
  main = main[which( !(main$V1 %in% c("-", ".")) ), ]
  N_pos = col_end-col_start+1-nrow(main_gaps)+cys_pre_cdr3_length
  identity = (N_pos-nrow(main))/N_pos
  if (identity==1){
    muts_main_nt_var = data.table(Mutations="",N = joined_filtered[which.max(joined_filtered$N),]$N)
    same_muts_as_main = joined_filtered
  }
  else{
    muts_main_nt_var = data.table(Mutations = paste(paste0(unlist(germline[,main$pos, with = FALSE]),                                               main$pos,
                                               main$V1), collapse = " "),
                      N = joined_filtered[which.max(joined_filtered$N), ]$N)
    
    who = mapply(function(x, y) { return(x == y) },
                 joined_filtered[, main$pos, with = FALSE],
                 main$V1)
    who = apply(who, 1, all)
    who = which(who)
    same_muts_as_main = joined_filtered[who, ]
  }
  #Extra mutations
  FR1toFR3 = same_muts_as_main[, -((col_end+1):ncol(same_muts_as_main))]
  to_count = FR1toFR3[, -(1:(col_start-1))] 
  to_count[to_count == "-"] = NA
  to_count[to_count == "."] = NA
  #to_count[to_count == "n"] = NA
  count_Nt = (rowSums(!is.na(to_count))) - nrow(main)
  same_muts_as_main$count_Nt = count_Nt
  not0 = same_muts_as_main[which(same_muts_as_main$count_Nt != 0), ]
  # TRY TO FIND LINKS IN THE SEQS WITH EXTRA MUTS --------------------------------------------------------------------
  # with dataset
  extralinks = not0[which(not0$N != 1), ]
  extralinks$Mutations = NA
  # get only the part of the alignment
  extralinks2 = extralinks[ ,col_start:col_end]
  # add the germline
  germline2 = germline[, col_start:col_end]
  # get the column names
  column_name = colnames(extralinks2)
  Mutations = mapply(function(x, y, pos) {
    return(paste0( y, pos, x))
  }, extralinks2, germline2, colnames(germline2))
  Mutations[str_detect(Mutations, "\\-")] = ""
  Mutations[str_detect(Mutations, "\\.")] = ""
  if (NCOL(Mutations)==1){
    Mutations = t(Mutations)
  }
  Mutations = apply(Mutations, 1, function(x) {
    return(paste(x, collapse = " "))
  })
  Mutations = str_squish(Mutations)
  extralinks$Mutations = Mutations
  if (identity<1){
    main_nt_variants = paste0(as.vector(unlist(germline[, main$pos, with = FALSE])),
                              main$pos,
                              main$V1)
    for(i in main_nt_variants) {
      extralinks$Mutations = str_remove_all(extralinks$Mutations, i)
    }
  }
  extralinks$Mutations = str_squish(extralinks$Mutations)
  extramutslinked = extralinks[, c("cluster_id", "N", "freq_cluster_id", "count_Nt", "Mutations"), with = FALSE]
  extramutslinked = extramutslinked[, .(cluster_id = paste(cluster_id, collapse = "+"),
                                        N = sum(N),
                                        count_Nt = unique(count_Nt)), by = Mutations]
  # eliminate cases with less than min_reads reads -----------------------------------------------------------------
  extramutslinked = extramutslinked[which(extramutslinked$N >= min_reads), ]
  # Error if no extra mutations
  if (nrow(extramutslinked)==0){
    mess = sprintf('The sample %s does not contain any extra mutations. The analysis has been stopped.',as.character(sample_id))
    f=file(paste0(save_path,'/warning_message.txt'))
    writeLines(mess, f)
    close(f)
    warning(mess)
    id_and_error_type = c(round(identity*100,2),TRUE,FALSE,FALSE)
    names(id_and_error_type) = c("identity","type_error_no_extra_nodes","type_error_max_length_1","type_error_no_connected_extra_nodes")
    save(id_and_error_type,file=paste0(save_path,'/id_and_error_type.Rda'))
    return(id_and_error_type)
  }
  # Error if maximal length of path is strictly smaller than 2
  if (max(extramutslinked$count_Nt)<2){
    mess = sprintf('The sample %s does not contain any pathways of length 2 or more. The analysis has been stopped.',as.character(sample_id))
    f=file(paste0(save_path,'/warning_message.txt'))
    writeLines(mess, f)
    close(f)
    warning(mess)
    id_and_error_type = c(round(identity*100,2),FALSE,TRUE,FALSE)
    names(id_and_error_type) = c("identity","type_error_no_extra_nodes","type_error_max_length_1","type_error_no_connected_extra_nodes")
    save(id_and_error_type,file=paste0(save_path,'/id_and_error_type.Rda'))
    return(id_and_error_type)
  }
  # EXTRA MUTATIONS -----------------------------------------------------------------  
  extra_connections = find_extra_connection(extramutslinked,include_jump=include_jump)
  if (nrow(extra_connections)==0){
    mess = sprintf('The sample %s does not contain any connected extra nodes. The analysis has been stopped.',as.character(sample_id))
    f=file(paste0(save_path,'/warning_message.txt'))
    writeLines(mess, f)
    close(f)
    warning(mess)
    id_and_error_type = c(round(identity*100,2),FALSE,FALSE,TRUE)
    names(id_and_error_type) = c("identity","type_error_no_extra_nodes","type_error_max_length_1","type_error_no_connected_extra_nodes")
    save(id_and_error_type,file=paste0(save_path,'/id_and_error_type.Rda'))
    return(id_and_error_type)
  }
  
  #get ordered unique values of numbers of extra mutations
  numbers2 = sort(unique(extramutslinked$count_Nt))
  
  #initiate all variables as empty lists
  extramut =list()
  extraseqs = list()

  k=1
  for(i in numbers2){
    filterextramuts = filter(extramutslinked, count_Nt == i)
    extramut[k] = filterextramuts %>% distinct(filterextramuts$Mutations) %>% nrow()  #-1 to not count the reference
    print(paste0(i, " extra muts = ", extramut[k]))
    extraseqs[k] = sum(filterextramuts$N, na.rm=TRUE) 
    print(paste0(i, " extra muts seqs = ", extraseqs[k]))
    k=k+1
    }
  
  extramutations = as.data.table(numbers2)
  extramutations$nt_var = as.character(extramut)
  extramutations$seqs = as.character(extraseqs)
  
  colnames(extramutations)[1]="#mutations"

  write.table(extramutations, paste0(save_path,"/extra_mutations_calculations_",sample_id,".txt"), sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # LESS MUTATIONS -----------------------------------------------------------------
  if (identity<1){
    less_expanded = get_less_mutations(joined_filtered,germline2,main,min_reads,col_start,col_end)
    if(nrow(less_expanded)== 0){
      less_expanded = c()
    }
    else{
    less_expanded2 = less_expanded[,2:3]
    less_expanded2 = aggregate(less_expanded2$N, by=list(less_expanded2$suppl_mut), FUN=sum)
    colnames(less_expanded2)[1]="#mutations"
    colnames(less_expanded2)[2]="N"
    write.table(less_expanded2, paste0(save_path,"/less_muts_calculations_",sample_id,".txt"), sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }}
  else{
    less_expanded = c()
  }
  main_node = data.table(Mutations = muts_main_nt_var$Mutations,
                         suppl_mut = 0,
                         N =muts_main_nt_var$N)
  evolution = rbind(less_expanded,main_node,extra_connections)
  evolution$level = "main"
  evolution[which(suppl_mut < 0), ]$level = "less"
  evolution[which(suppl_mut > 0), ]$level = "additional"
  evolution$Mutations = str_squish(evolution$Mutations)
  write.table(evolution,paste0(save_path,"/evolution-file-",sample_id,"_threshold",min_reads,".txt"), sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
  evolution = evolution[,c("Mutations", "suppl_mut", "level", "N"), with = FALSE]
  evolution = unique(evolution)
  # create adjacency matrix -----------------------------
  evolution = evolution[which(!is.na(evolution$Mutations)), ]
  mutation_pos = str_split(evolution$Mutations, " ")
  names(mutation_pos) = evolution$Mutations
  adjacency = lapply(mutation_pos, function(x) {
    out = lapply(mutation_pos, function(y, xt) {
      return( max(length(which( !(xt %in% y) )),
                  length(which( !(y %in% xt) )) ) )
    }, x)
    out = data.table(from = paste(x, collapse = " "),
                     to = names(out),
                     diff_n = unlist(out))
    return(out)
  })
  adjacency = rbindlist(adjacency)
  who = base::match(adjacency$from, evolution$Mutations)
  adjacency$suppl_mu_f = evolution[who, ]$suppl_mut
  adjacency$length_f = unlist(lapply(str_split(adjacency$from, " "), length))
  adjacency$n_reads_f = evolution[who, ]$N
  who = base::match(adjacency$to, evolution$Mutations)
  adjacency$suppl_mu_t = evolution[who, ]$suppl_mut
  adjacency$length_t = unlist(lapply(str_split(adjacency$to, " "), length))
  adjacency$n_reads_t = evolution[who, ]$N
  # keep evolution connections -------------------
  adjacency = adjacency[which( adjacency$suppl_mu_f < adjacency$suppl_mu_t ), ]
  #who = which(adjacency$suppl_mu_f == 0 & adjacency$suppl_mu_t > 0)
  #adjacency[who, ]$diff_n = adjacency[who, ]$suppl_mu_t
  if (include_jump){
    adjacency = adjacency[which(adjacency$suppl_mu_f==0 | adjacency$diff_n == (adjacency$length_t - adjacency$length_f) | (adjacency$suppl_mu_f<0 & adjacency$suppl_mu_t <= 0)), ]
    ind2keep = logical(nrow(adjacency))
    checked_ind = list()
    for (i in 1:nrow(adjacency)){
      if (!(i %in% checked_ind)){
        who = which(adjacency$to==adjacency$to[i])
        checked_ind = append(checked_ind,who)
        ind2keep[who[adjacency$diff_n[who]==min(adjacency$diff_n[who])]]=TRUE
      }
    }
    adjacency = adjacency[ind2keep,]
  }
  else{
    adjacency = adjacency[which(adjacency$diff_n <= 1), ]
    adjacency = adjacency[which(abs(adjacency$suppl_mu_f - adjacency$suppl_mu_t) <= 1), ]
  }
  # create connection matrix --------------------
  #edges = adjacency[which((adjacency$diff_n != 0) | (adjacency$from == adjacency$to)), ]
  edges = adjacency[which(adjacency$suppl_mu_f <= adjacency$suppl_mu_t), ]
  edges = edges[order(suppl_mu_f), ]
  edges = edges[which(edges$from != edges$to), ]
  who = which(!(edges$to %in% edges$from) & edges$suppl_mu_t < 0)
  if(length(who) > 0) {
    edges = edges[-who, ]
  }
  nodes = rbind(edges[, c("from", "n_reads_f", "suppl_mu_f"), with = FALSE], 
                edges[, c("to", "n_reads_t", "suppl_mu_t"), with = FALSE], use.names = FALSE)
  nodes = unique(nodes)
  colnames(nodes) = c("Mutations", "N", "suppl_mu")
  nodes$id = paste0("N", 1:nrow(nodes))
  nodes = nodes[, c("id", "Mutations", "N", "suppl_mu"), with = FALSE]
  who = match(nodes$Mutations, evolution$Mutations)
  nodes$level = evolution[who, ]$level
  who = match(edges$from, nodes$Mutations)
  edges$from = nodes[who,]$id
  who = match(edges$to, nodes$Mutations)
  edges$to = nodes[who,]$id
  
  edges = edges[, c("from", "to"), with = FALSE]
  graph_info = data.table(
    sample_id=sample_id,
    main_nt_var_identity=round(identity*100,2),
    convergence_score= double(1),
    nb_reads_most_relevant_pathway=integer(1),
    nb_reads_main_nt_var=integer(1),
    most_relevant_pathway_score=double(1),
    nb_nodes_most_relevant_pathway=integer(1),
    max_path_length=integer(1),
    max_muts_length=integer(1),
    end_nodes_density=double(1),
    nb_end_nodes=integer(1),
    nb_extra_nodes=integer(1),
    nb_reads_tot=integer(1),
    avg_degree=double(1),
    avg_distance=integer(1),
    stringsAsFactors=FALSE)

  block_table = get_main_block(nodes,edges)
  graph_info$most_relevant_pathway_score = as.double(block_table$most_relevant_pathway_score[1])
  graph_info$nb_nodes_most_relevant_pathway = as.integer(block_table$nb_nodes_most_relevant_pathway[1])
  graph_info$nb_reads_most_relevant_pathway = as.integer(block_table$nb_reads_most_relevant_pathway[1])
  graph_info$nb_reads_main_nt_var = as.integer(block_table$nb_reads_main_nt_var[1])
  graph_info$convergence_score = as.double(block_table$convergence_score[1])
  graph_info$nb_reads_tot = as.integer(sum(nodes$N))
  graph_info$max_path_length = max(block_table$path_length)
  write.table(block_table, paste0(save_path,'/block_table_',sample_id,'.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  # nodes = get_nodes_aa_muts(nodes,germline2)
  nodes = get_nodes_aa_muts(nodes,germline[,col_pos1:col_end])
  aa_muts_weight = get_aa_muts_weight(nodes)
  aa_muts_weight_main_variant = get_aa_muts_weight(nodes,with_main_nt_variant=TRUE)
  #nodesO = nodes
  #edgesO = edges
  if (nrow(aa_muts_weight)>0){
    write.table(aa_muts_weight, paste0(save_path,'/aa_muts_weight_',sample_id,'.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  }
  if (nrow(aa_muts_weight_main_variant)>0){
    write.table(aa_muts_weight_main_variant, paste0(save_path,'/aa_muts_weight_main_variant_',sample_id,'.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  }
  save(nodes,file=paste0(save_path,'/nodes_all_',sample_id,'.Rda'))
  save(edges,file=paste0(save_path,'/edges_all_',sample_id,'.Rda'))
  edges = select_main_block(block_table,nodes,edges,ind2take=1,include_max_extra_muts=TRUE,include_max_path_length=TRUE)
  nodes = as.data.frame(edges[1])
  edges = as.data.frame(edges[2])
  save(nodes,file=paste0(save_path,'/nodes_graph_',sample_id,'.Rda'))
  save(edges,file=paste0(save_path,'/edges_graph_',sample_id,'.Rda'))
  write.table(nodes, paste0(save_path,'/nodes_graph_',sample_id,'.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  
  nodes$Mutations_aa[nodes$suppl_mu==0]=""
  tidy_net = tbl_graph(nodes = nodes, 
                       edges = edges, 
                       directed = TRUE)
  
  set.seed(10)
  color_start = 3
  max_level1 = RColorBrewer::brewer.pal.info["Reds",]$maxcolors
  max_level2 = RColorBrewer::brewer.pal.info["Greens",]$maxcolors
  min_level1 = RColorBrewer::brewer.pal.info["Blues",]$maxcolors
  min_level2 = RColorBrewer::brewer.pal.info["Purples",]$maxcolors
  colors_plus1 = RColorBrewer::brewer.pal(max_level1, "Reds")[color_start:max_level1]
  colors_plus2 = RColorBrewer::brewer.pal(max_level2, "Greens")[color_start:max_level2]
  colors_minus1 = RColorBrewer::brewer.pal(min_level1, "Blues")[color_start:min_level1]
  colors_minus2 = RColorBrewer::brewer.pal(min_level2, "Purples")[color_start:min_level2]
  colors_plus = vector(class(colors_plus1), length(c(colors_plus1, colors_plus2)))
  colors_plus[c(TRUE, FALSE)] = colors_plus1
  colors_plus[c(FALSE, TRUE)] = colors_plus2
  colors_minus = vector(class(colors_minus1), length(c(colors_minus1, colors_minus2)))
  colors_minus[c(TRUE, FALSE)] = colors_minus1
  colors_minus[c(FALSE, TRUE)] = colors_minus2
  color0 = '#005F6A'
  color_values = append(append(rev(colors_minus),color0),colors_plus)
  max_level = max_level1+max_level2-(2*color_start-2)
  min_level = min_level1+min_level2-(2*color_start-2)
  names(color_values)=as.factor((-min_level):max_level)
  color_default = '#8B8989'
  if (min(nodes$suppl_mu)<(-min_level)){
    n2add = abs(min(nodes$suppl_mu))-min_level
    color2add = rep(color_default,n2add)
    names(color2add) = as.factor(min(nodes$suppl_mu):(-min_level-1))
    color_values = append(color2add,color_values)
    warning(sprintf("In sample %s, there are more negative mutations (%i) than colors for the nodes (%i). A default color is used for those nodes.",as.character(sample_id),abs(min(nodes$suppl_mu)),min_level))
  }
  if (max(nodes$suppl_mu)>max_level){
    n2add = max(nodes$suppl_mu)-max_level
    color2add = rep(color_default,n2add)
    names(color2add) = as.factor((max_level+1):max(nodes$suppl_mu))
    color_values = append(color_values,color2add)
    warning(sprintf("In sample %s, there are more positive mutations (%i) than colors for the nodes (%i). A default color is used for those nodes.",as.character(sample_id),max(nodes$suppl_mu),max_level))
  }
  
  options(ggrepel.max.overlaps = Inf)
  rr = which(nodes$suppl_mu==min(nodes$suppl_mu))
  if (nodes_size_scaling) {
    if (include_aa_muts) {
      ggraph(tidy_net, layout = layout_as_tree(tidy_net,root=rr,mode='all')) +
        geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                      start_cap = circle(2, 'mm'),
                      end_cap = circle(3, 'mm'),
                      alpha = 0.1) + 
        geom_node_point(aes(size = log10(N),
                          color = as.factor(suppl_mu)),
                      alpha = 0.7) +
        geom_node_text(aes(label = Mutations_aa), repel=TRUE,size=2,check_overlap=FALSE) +
        scale_color_manual(name = "Extra mutations", values=color_values,drop=TRUE,limits=force) +
        theme(panel.background = element_blank(),
              legend.background = element_blank()) +
        labs(subtitle = paste("Graph Network", sample_id))
    }
    else {
      ggraph(tidy_net, layout = layout_as_tree(tidy_net,root=rr,mode='all')) +
        geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                       start_cap = circle(2, 'mm'),
                       end_cap = circle(3, 'mm'),
                       alpha = 0.1) + 
        geom_node_point(aes(size = log10(N),
                            color = as.factor(suppl_mu)),
                        alpha = 0.7) +
        scale_color_manual(name = "Extra mutations", values=color_values,drop=TRUE,limits=force) +
        theme(panel.background = element_blank(),
              legend.background = element_blank()) +
        labs(subtitle = paste("Graph Network", sample_id))
    }
  }
  else {
    if (include_aa_muts) {
      ggraph(tidy_net, layout = layout_as_tree(tidy_net,root=rr,mode='all')) +
        geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                      start_cap = circle(2, 'mm'),
                      end_cap = circle(3, 'mm'),
                      alpha = 0.1) + 
        geom_node_point(aes(color = as.factor(suppl_mu)),
                        alpha = 0.7) +
        geom_node_text(aes(label = Mutations_aa), repel=TRUE,size=2,check_overlap=FALSE) +
        scale_color_manual(name = "Extra mutations", values=color_values,drop=TRUE,limits=force) +
        theme(panel.background = element_blank(),
              legend.background = element_blank()) +
        labs(subtitle = paste("Graph Network", sample_id))
    }
    else {
      ggraph(tidy_net, layout = layout_as_tree(tidy_net,root=rr,mode='all')) +
        geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                       start_cap = circle(2, 'mm'),
                       end_cap = circle(3, 'mm'),
                       alpha = 0.1) + 
        geom_node_point(aes(color = as.factor(suppl_mu)),
                        alpha = 0.7) +
        scale_color_manual(name = "Extra mutations", values=color_values,drop=TRUE,limits=force) +
        theme(panel.background = element_blank(),
              legend.background = element_blank()) +
        labs(subtitle = paste("Graph Network", sample_id))
    }
  }
  ggsave(paste0(save_path,'/',sample_id,'_ID.pdf'), width = 9.0, height = 9.0, units="in")
  #start graph metrics calculation -----------------------------
  d = degree(tidy_net, mode = "all", loops = FALSE, normalized = FALSE)
  avg_degree = mean(d)
  avg_dist = mean_distance(tidy_net, 
                           directed = TRUE,
                           unconnected = TRUE)
  tidy_net2 = as.data.frame(tidy_net)
  extra = tidy_net2[which(tidy_net2$suppl_mu > 0),]
  nb_extra_nodes = nrow(extra)
  nb_end_nodes = length(intersect(extra$id,setdiff(edges$to,edges$from)))
  end_nodes_density = nb_end_nodes / nb_extra_nodes
  max_muts_length = max(tidy_net2$suppl_mu)
  graph_info$avg_degree = round(avg_degree,3)
  graph_info$avg_distance = as.integer(avg_dist)
  graph_info$nb_extra_nodes = as.integer(nb_extra_nodes)
  graph_info$max_muts_length = as.integer(max_muts_length)
  graph_info$nb_end_nodes = as.integer(nb_end_nodes)
  graph_info$end_nodes_density = round(end_nodes_density,3)
  write.table(graph_info, paste0(save_path,'/graph_info_',sample_id,'.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  return(graph_info)
}

find_extra_connection = function(extramutslinked,include_jump=TRUE){
  extra_connections = list()
  max_extra_muts = max(extramutslinked$count_Nt)
  if (include_jump){
    max_jump = max_extra_muts
    
  }
  else{
    max_jump = 1
  }
  k = 1
  for (jump_i in 1:max_jump){
    for(i in 1:nrow(extramutslinked)) {
      x = str_split(extramutslinked[i, ]$Mutations, "\\ ")
      x = unlist(x)
      if ((length(x) + jump_i)<=max_extra_muts){
        temp = extramutslinked[which( extramutslinked$count_Nt == (length(x) + jump_i) ), ]
      }
      else{
        next
      }
      if(nrow(temp) == 0) {
        next
      }
      Mutations = str_split(temp$Mutations, "\\ ")
      who = lapply(Mutations, function(y) {
        return(all(x %in% y))
        })
      who = unlist(who)
      who = which(who)
      if(length(who) == 0) {
        next
      }
      temp = temp[who, ]
      extra_connections[[k]] = rbind(data.table(Mutations = extramutslinked[i, ]$Mutations,suppl_mut = extramutslinked[i, ]$count_Nt,N = extramutslinked[i, ]$N),
                                     data.table(Mutations = temp$Mutations,suppl_mut = temp$count_Nt,N = temp$N))
      k = k + 1
    }
  }
  extra_connections = rbindlist(extra_connections)
  extra_connections = unique(extra_connections)
  if (nrow(extra_connections)>0){
    extra_connections = extra_connections[order(suppl_mut,-N,Mutations),]
  }
  return(extra_connections)
}

get_less_mutations = function(joined_filtered,germline2,main,min_reads=10,col_start=5,col_end=313){
  # joined filtered no CDR3
  jfnc = joined_filtered[,-((col_end+1):ncol(joined_filtered))]
  who = colnames(jfnc)[which(!(colnames(jfnc) %in% main$pos))]
  jfnc_not_main = jfnc[, who, with = FALSE]
  jfnc_not_main = jfnc_not_main[, col_start:ncol(jfnc_not_main)]
  jfnc_not_main[jfnc_not_main == '-'] = NA
  jfnc_not_main[jfnc_not_main == '.'] = NA
  who = apply(jfnc_not_main, 1, function(x){ 
    return(all(is.na(x)))
  })
  jfnc_not_main = jfnc_not_main[who, ]
  jfnc = jfnc[who, ]
  jfnc_main = jfnc[, main$pos, with = FALSE]
  jfnc_main[jfnc_main == '-'] = NA
  who = apply(jfnc_main, 1, function(x){ 
    return(any(is.na(x)))
  })
  jfnc_main = jfnc_main[which(who), ]
  jfnc = jfnc[which(who), ]
  less_expanded = jfnc[which(jfnc$N > 1), ]
  Mutations = mapply(function(x, y, pos) {
    
    return(paste0( y, pos, x))
    
  }, less_expanded[, col_start:ncol(less_expanded)], germline2, colnames(germline2))
  Mutations[str_detect(Mutations, "\\.")] = ""
  if (NCOL(Mutations)==1){
    Mutations = t(Mutations)
  }
  Mutations = Mutations[, main$pos, drop = FALSE]
  Mutations = apply(Mutations, 1, function(x) {
    return(paste(x, collapse = " "))
  })
  Mutations = str_squish(Mutations)
  less_expanded$Mutations = Mutations
  less_expanded = less_expanded[, .(cluster_id = paste(cluster_id, collapse = "+"),
                                    N = sum(N)), by = Mutations]
  less_expanded = less_expanded[which(less_expanded$N >= min_reads), ]
  less_expanded$suppl_mut = -1 * str_count(less_expanded$Mutations, "\\-")
  less_expanded = as.data.table(less_expanded)
  less_expanded = less_expanded[ , c("Mutations", "suppl_mut", "N"), with = FALSE]
  less_expanded = less_expanded[order(suppl_mut,-N,Mutations),]
  return(less_expanded)
}


