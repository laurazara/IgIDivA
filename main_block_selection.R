recursive_search = function(node_curr,edges,main_node,nodes_in_block){
  if (node_curr!=main_node){
    next_nodes = edges$from[edges$to==node_curr]
    for (nn in next_nodes){
      nodes_in_block = append(nodes_in_block,nn)
    }
    for (nn in next_nodes){
      nodes_in_block = append(nodes_in_block,recursive_search(nn,edges,main_node,nodes_in_block))
    }
    return(unique(nodes_in_block))
  }
  else {
    return(unique(nodes_in_block))
  }
}

get_main_block = function(nodes,edges){
  end_node = setdiff(edges$to,edges$from)
  end_node = setdiff(end_node,nodes$id[nodes$suppl_mu<1])
  num_block = length(end_node)
  main_node = nodes$id[nodes$suppl_mu==0]
  block_table = data.frame(
                   end_node=end_node,
                   most_relevant_pathway_score=double(num_block),
                   convergence_score=double(num_block),
                   nb_reads_most_relevant_pathway=integer(num_block),
                   nb_reads_main_nt_var=integer(num_block),
                   nb_nodes_most_relevant_pathway=integer(num_block),
                   nodes=character(num_block),
                   path_length=integer(num_block),
                   muts_length=integer(num_block),
                   stringsAsFactors=FALSE)
  total_weight = sum(nodes$N[nodes$suppl_mu>0])
  main_nodes_reads = nodes$N[nodes$id==main_node]
  for (i in 1:num_block){
    node_curr = end_node[i]
    nodes_in_block = recursive_search(node_curr,edges,main_node,character())
    nodes_in_block = nodes_in_block[nodes_in_block!=main_node]
    w = nodes$N[nodes$id==node_curr]
    l = nodes$suppl_mu[nodes$id==node_curr]
    for (n in nodes_in_block){
      w = w + nodes$N[nodes$id==n]
      l = append(l,nodes$suppl_mu[nodes$id==n])
    }
    block_table$most_relevant_pathway_score[i] = round(as.double(w)/total_weight,digits=3)
    block_table$convergence_score[i] = round(w/main_nodes_reads,3)
    block_table$nb_reads_most_relevant_pathway[i] = as.integer(w)
    block_table$nb_reads_main_nt_var[i] = main_nodes_reads
    block_table$nb_nodes_most_relevant_pathway[i] = length(nodes_in_block)
    block_table$nodes[i] = paste(nodes_in_block,collapse=',')
    block_table$path_length[i] = length(unique(l))
    block_table$muts_length[i] = nodes$suppl_mu[nodes$id==node_curr]
  }
  block_table = block_table[
    with(block_table, order(most_relevant_pathway_score,nb_nodes_most_relevant_pathway,decreasing = TRUE)),
    ]
  return(block_table)
}

select_main_block = function(block_table,nodes,edges,ind2take=0,n_min=1,min_score=0.1,include_max_extra_muts=TRUE,include_max_path_length=TRUE){
  if (ind2take==0){
    ind2take = block_table$most_relevant_pathway_score>min_score
    n_min = min(nrow(block_table),n_min)
    ind2take[1:n_min] = TRUE
    ind2take = which(ind2take)
  }
  if (include_max_extra_muts){
    ind2take = unique(append(ind2take,which(block_table$muts_length==max(block_table$muts_length))))
  }
  if (include_max_path_length){
    ind2take = unique(append(ind2take,which(block_table$path_length==max(block_table$path_length))))
  }
  unique_nodes = strsplit(block_table$nodes[ind2take],',')
  unique_nodes = unique(do.call(c,unique_nodes))
  unique_nodes = append(unique_nodes,nodes$id[nodes$suppl_mu<=0])
  unique_nodes = append(unique_nodes,block_table$end_node[ind2take])
  nodes2take = logical(length = nrow(nodes))
  edges2take_f = logical(length = nrow(edges))
  edges2take_t = edges2take_f
  for (n in unique_nodes){
    nodes2take = nodes2take | nodes$id==n
    edges2take_f = edges2take_f | edges$from==n
    edges2take_t = edges2take_t | edges$to==n
  }
  nodes = nodes[nodes2take,]
  edges = edges[edges2take_f & edges2take_t,]
  return(list(nodes,edges))
}