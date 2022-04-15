source("id_network_analysis.R")

iterate_do_graph = function(argument_file,save_path = getwd(),include_jump=TRUE,col_start = 5,col_end = 313,min_reads = 10,highly_sim_clonos = c(1),nodes_size_scaling = TRUE,include_aa_muts = TRUE){
  output_path = paste0(save_path, "/Output")
  dir.create(output_path,showWarnings = FALSE)
  args = data.table::fread(argument_file,
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE)
  
  
  n = nrow(args)
  # no_muts_df <- data.frame(sample_id=character(n))
  aa_muts_all = data.table()
  aa_muts_main_variant_all = data.table()
  metric_table = data.frame(
    sample_id=character(n),
    main_nt_var_identity=double(n),
    convergence_score= double(n),
    nb_reads_most_relevant_pathway=integer(n),
    nb_reads_main_nt_var=integer(n),
    most_relevant_pathway_score=double(n),
    nb_nodes_most_relevant_pathway=integer(n),
    max_path_length=integer(n),
    max_muts_length=integer(n),
    end_nodes_density=double(n),
    nb_end_nodes=integer(n),
    nb_extra_nodes=integer(n),
    nb_reads_tot=integer(n),
    avg_degree=double(n),
    avg_distance=integer(n),
    error_type=character(n),
    stringsAsFactors=FALSE)
  n_col = ncol(metric_table)
  for (i in 1:n){
    out_file_name = paste0(save_path,'/',args$sample_id[i],'/graph_info_',args$sample_id[i],'.txt')
    metric_table[i,1] = args$sample_id[i]
    if (file.exists(out_file_name)){
      graph_info = data.table::fread(out_file_name,
                                     header = TRUE,
                                     sep = "\t",
                                     stringsAsFactors = FALSE)
      metric_table[i,1:(n_col-1)] = graph_info
    }
    else {
      error_file = paste0(save_path,'/',args$sample_id[i],'/id_and_error_type.Rda')
      if (file.exists(error_file)){
        load(error_file)
        graph_info = id_and_error_type
        # temp <- args$sample_id[i]
        # temp <- as.data.frame(temp)
        # colnames(temp) <- "sample_id"
        # no_muts_df <- rbind(no_muts_df, temp)
      }
      else if (!file.exists(paste0(save_path,'/',args$sample_id[i]))){
        graph_info = doGraph(args$highly_sim_clonos_file[i],args$grouped_alignment_file[i],args$sample_id[i],save_path=save_path,include_jump=include_jump,col_start = col_start, col_end = col_end, min_reads = min_reads, highly_sim_clonos = highly_sim_clonos, nodes_size_scaling = nodes_size_scaling, include_aa_muts = include_aa_muts)
      }
      else {
        graph_info = c(0/0,FALSE,FALSE,FALSE)
      }
      if (typeof(graph_info[1])=="double"){
        metric_table[i,3:(n_col-1)] = 0/0
        metric_table[i,2] = graph_info[1]
        type_error = as.logical(graph_info[2:length(graph_info)])
        names(type_error) = names(graph_info[2:length(graph_info)])
        if (type_error[1]){
          metric_table$max_muts_length[i] = 0
        }
        else if (type_error[2]){
          metric_table$max_muts_length[i] = 1
        }
        if (any(type_error)){
          metric_table$error_type[i] = names(type_error)[which(type_error)]
        }
        else{
          metric_table$error_type[i] = "sample_output_folder_exits_without_result_or_error_files"
        }
      }
      else{
        metric_table[i,1:(n_col-1)] = graph_info
      }
    }
    file_aa_muts = paste0(save_path,'/',args$sample_id[i],'/aa_muts_weight_',args$sample_id[i],'.txt')
    file_aa_muts_main_variant = paste0(save_path,'/',args$sample_id[i],'/aa_muts_weight_main_variant_',args$sample_id[i],'.txt')
    if (file.exists(file_aa_muts)){
      aa_muts = data.table::fread(file_aa_muts,
                               header = TRUE,
                               sep = "\t",
                               stringsAsFactors = FALSE)
      aa_muts$id = rep(args$sample_id[i],nrow(aa_muts))
      aa_muts_all = rbind(aa_muts_all,aa_muts)
    }
    if (file.exists(file_aa_muts_main_variant)){
      aa_muts_main_variant = data.table::fread(file_aa_muts_main_variant,
                                  header = TRUE,
                                  sep = "\t",
                                  stringsAsFactors = FALSE)
      aa_muts_main_variant$id = rep(args$sample_id[i],nrow(aa_muts_main_variant))
      aa_muts_main_variant_all = rbind(aa_muts_main_variant_all,aa_muts_main_variant)
    }
    
  }
  write.table(metric_table, paste0(output_path,'/metric_table_all.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  write.table(aa_muts_all, paste0(output_path,'/aa_muts_weight.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  write.table(aa_muts_main_variant_all, paste0(output_path,'/aa_muts_weight_main_variant.txt'), sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote = FALSE)
  
  
  discarded_samples <- data.table::fread(paste0(output_path,'/metric_table_all.txt'), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  discarded_samples <- na.omit(discarded_samples, invert=TRUE)
  discarded_samples <- as.data.frame(discarded_samples$sample_id)
  colnames(discarded_samples) <- "sample_id"
  write.table(discarded_samples, paste0(output_path,"/discarded_samples_table.txt"), sep = "\t", dec = ".", 
              row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)
  
  return(metric_table)
}
