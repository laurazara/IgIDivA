# Define functions
create_summary_calculations_table <- function(sample, save_path) {
    summary_calculations <- data.table::fread(paste0(save_path, "/", sample, "/summary-calculations_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(summary_calculations)
}

create_extra_mutations_calculations_table <- function(sample, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/extra_mutations_calculations_", sample, ".txt"))){
        extra_mutations_calculations <- data.table::fread(paste0(save_path, "/", sample, "/extra_mutations_calculations_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        return(extra_mutations_calculations)
    }else{
        return()
    }
}

create_less_mutations_calculations_table <- function(sample, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/less_muts_calculations_", sample, ".txt"))){
        less_mutations_calculations <- data.table::fread(paste0(save_path, "/", sample, "/less_muts_calculations_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        return(less_mutations_calculations)
    }else{
        return()
    }
}

create_mutations_table <- function(sample, min_reads, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/evolution-file-", sample, "_threshold", min_reads, ".txt"))){
        mutations <- data.table::fread(paste0(save_path, "/", sample, "/evolution-file-", sample, "_threshold", min_reads, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        return(mutations)
    }else{
        return()
    }
}

create_aa_mutations_main_var_table <- function(sample, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/aa_muts_weight_main_variant_", sample, ".txt"))){
        aa_mutations_main_var <- data.table::fread(paste0(save_path, "/", sample, "/aa_muts_weight_main_variant_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        return(aa_mutations_main_var)
    }else{
        return()
    }
}

create_aa_mutations_table <- function(sample, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/aa_muts_weight_", sample, ".txt"))){
        aa_mutations <- data.table::fread(paste0(save_path, "/", sample, "/aa_muts_weight_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        return(aa_mutations)
    }else{
        return()
    }
}

create_graph_metrics_table <- function(sample, include_metrics, save_path) {
    
    if(file.exists(paste0(save_path, "/", sample, "/graph_info_", sample, ".txt"))){
        graph_metrics <- data.table::fread(paste0(save_path, "/", sample, "/graph_info_", sample, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        metrics <- unlist(strsplit(include_metrics, ","))
        
        if (!("Main variant identity" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-main_nt_var_identity)
        }
        if (!("Relative convergence (reads)" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-convergence_score, -nb_reads_most_relevant_pathway, -nb_reads_main_nt_var)
        }
        if (!("Most relevant pathway score" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-most_relevant_pathway_score)
        }
        if (!("Most relevant pathway (nodes)" %in%  metrics)) {
            graph_metrics <- graph_metrics %>% select(-nb_nodes_most_relevant_pathway)
        }
        if (!("End nodes density" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-end_nodes_density, -nb_end_nodes, -nb_extra_nodes)
        }
        if (!("Max path length" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-max_path_length)
        }
        if (!("Max mutations path length" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-max_muts_length)
        }
        if (!("Total reads" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-nb_reads_tot)
        }
        if (!("Average degree" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-avg_degree)
        }
        if (!("Average distance" %in% metrics)) {
            graph_metrics <- graph_metrics %>% select(-avg_distance)
        }
        return(graph_metrics)
        
    }else{
        return()
    }
    
}

create_graph_network <- function(sample, save_path) {
    if(file.exists(paste0(save_path, "/", sample, "/", sample, "_ID.pdf"))){
        pdf_convert(paste0(save_path, "/", sample, "/", sample, "_ID.pdf"), format = "png", filenames = paste0(save_path, "/", sample, "/", sample, "_ID.png"), dpi = 300)
        list(src = paste0(save_path, "/", sample, "/", sample, "_ID.png"), width = 800, height = 800)
    }else{
        return()
    }
}

choose_comparison_metrics <- function(compare_metrics) {
    comparison_metrics <- c("convergence_score", "end_nodes_density", "max_path_length", "max_muts_length", "avg_degree", "avg_distance")
    metrics <- unlist(strsplit(compare_metrics, ","))
    if (!("Convergence score" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "convergence_score"]
    }
    if (!("End nodes density" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "end_nodes_density"]
    }
    if (!("Max path length" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "max_path_length"]
    }
    if (!("Max mutations length" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "max_muts_length"]
    }
    if (!("Average degree" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "avg_degree"]
    }
    if (!("Average distance" %in% metrics)) {
        comparison_metrics <- comparison_metrics[!comparison_metrics %in% "avg_distance"]
    }
    return(comparison_metrics)
}
