library(data.table)
library(stringr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(rstatix)

perform_statistical_analysis <- function(
    groups_file, 
    final_metric_table, 
    comparison_metrics = c(
        "convergence_score", 
        "end_nodes_density", 
        "max_path_length", 
        "max_muts_length", 
        "avg_degree", 
        "avg_distance"
    ), 
    compare_by = "group_name",
    p_threshold = 0.05,
    adjust = "no"
) {
  
  dir.create(paste0(getwd(), "/Comparisons"), showWarnings = FALSE)
    
    if (adjust=='no'){
        adj_text = ''
    }
    else{
        adj_text = 'adjusted'
    }
    
  
  # read inputs -------------------------------------
    
    if(file.exists("groups_file")){
        sample_metadata = fread(
            groups_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE
        )
    }
  

  
  data_metrics = fread(
    final_metric_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
  
  metrics_to_compare = comparison_metrics
  
  groups_to_compare = compare_by
  
  # merge data and sample metadata ------------------------------
  
  metrics_to_compare = metrics_to_compare[which(
    metrics_to_compare %in% colnames(data_metrics)
  )]
  
  data_metrics = merge(data_metrics, sample_metadata, by = "sample_id")

  
  # kruskal statistics ------------------------------------
  
  krusk = list()
  
  
  for(g in groups_to_compare) {
    
    
    for(m in metrics_to_compare) {
      
      tmp = data_metrics[which(
        !is.na(data_metrics[[m]]) & !is.na(data_metrics[[g]])
      ), ]
      
      
      tmp = kruskal_test(
        data = tmp,
        formula = as.formula( paste0(m, " ~ ", g) )
      )
      
      tmp = data.table(
        "comparison" = paste0(m, " ~ ", g),
        "p.val" = tmp$p
      )
      
      
      
      krusk[[ paste0(m, " ~ ", g) ]] = tmp
      
    }
    
    
  }
  
  krusk = rbindlist(krusk)
  
  # plot results ----------------------------
  
  message("Computing pairwise stats")
  
  plots_list = list()
  
  for(i in groups_to_compare) {
    
    for(j in metrics_to_compare) {
      
      tmp = data_metrics[which(
        !is.na(data_metrics[[ j ]]) & !is.na(data_metrics[[ i ]])
      ), ]
      
      tmp = tmp[, c(i, j), with = FALSE]
      
      colnames(tmp) = c("group", "metric")
      
      if(length(unique(tmp$metric)) <=  1) {
          
          next
          
      }
      
      # message(c("\nAdjust value = ", adjust, "\n", "Group: ", i, " Metric: ", j))
      
      if( adjust == "no" | length(unique(tmp$group)) <= 2 ){
          
          combinations = wilcox_test(
              data = tmp,
              formula = metric ~ group,
              p.adjust.method = "none"
          )
          
          combinations = setDT(combinations)
          
          combinations = combinations[which(combinations$p <= p_threshold), ]
          
          combinations$p.value = combinations$p
          
      } else {
          
          combinations = wilcox_test(
              data = tmp,
              formula = metric ~ group,
              p.adjust.method = "holm"
          )
          
          combinations = setDT(combinations)
          
          combinations = combinations[which(combinations$p.adj <= p_threshold), ]
          
          combinations$p.value = combinations$p.adj
      }
      
      if(nrow(combinations) == 0){
          
        s2 = tmp[, by = group, .N]
        
        tmp = merge(tmp, s2, by = "group")
        
        tmp$xlabel = paste0(
            tmp$group, "\n",
            "(n = ", tmp$N, ")"
        )

        gr = ggplot(data = tmp,
                    aes(x = xlabel, y = metric)) +


          geom_jitter(aes(col = group),
                      shape = 1,
                      width = 0.2,
                      height = 0) +

          geom_boxplot(aes(fill = group),
                       alpha = 0.5,
                       width = 0.5,
                       outlier.shape = NA) +

          stat_summary(fun = mean,
                       col = 'red',
                       geom = 'point',
                       shape = 20,
                       size = 8) +

          scale_fill_nejm() +
          scale_color_nejm() +

          theme_minimal() +

          theme(
              legend.position = "none",
              
              axis.title.x = element_blank(),
              
              panel.grid = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_line()
          ) +

          labs(
              y = j,
              # x = "Groups",
              caption = paste0(
                  "Kruskal-Wallis, p-val = ",
                  format(
                      krusk[which(
                          krusk$comparison == paste0(j, " ~ ", i)
                      ), ]$p.val, 
                      scientific = TRUE, 
                      digits = 4 
                  )
              )
          ) +
          theme(text = element_text(size = 18))

        ggsave(
          filename = paste0(getwd(), "/Comparisons/", j, "_vs_", i, ".pdf"),
          plot = gr, width = 9.0, height = 9.0, units = "in"
        )

        plots_list[[ paste0(j, "_vs_", i) ]] = gr

      } else {

        combinations = add_y_position(combinations)
        
        s2 = tmp[, by = group, .N]
        
        tmp = merge(tmp, s2, by = "group")
        
        tmp$xlabel = paste0(
            tmp$group, "\n",
            "(n = ", tmp$N, ")"
        )
        
        who = match(combinations$group1, tmp$group)
        combinations$group1 = tmp[who, ]$xlabel
        
        who = match(combinations$group2, tmp$group)
        combinations$group2 = tmp[who, ]$xlabel

        gr = ggplot(data = tmp,
                    aes(x = xlabel, y = metric)) +


          geom_jitter(aes(col = group),
                      shape = 1,
                      width = 0.2,
                      height = 0) +

          geom_boxplot(aes(fill = group),
                       alpha = 0.5,
                       width = 0.5,
                       outlier.shape = NA) +

          stat_summary(fun = mean,
                       col = 'red',
                       geom = 'point',
                       shape = 20,
                       size = 8) +

          stat_pvalue_manual(combinations, label = "p = {p.value}",
                             tip.length = 0.01) +

          scale_fill_nejm() +
          scale_color_nejm() +

          theme_minimal() +

          theme(
              legend.position = "none",
              
              axis.title.x = element_blank(),
              
              panel.grid = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_line()
          ) +

          labs(
              y = j,
              # x = "Groups",
              caption = paste0(
                'p:',adj_text,' Wilcoxon Rank Sum p-value\n',
                "Kruskal-Wallis, p-val = ",
                  format(
                      krusk[which(
                          krusk$comparison == paste0(j, " ~ ", i)
                      ), ]$p.val, 
                      scientific = TRUE, 
                      digits=4 
                  )
              )
          ) +
          theme(text = element_text(size = 18))

        ggsave(
          filename = paste0(getwd(), "/Comparisons/", j, "_vs_", i, ".pdf"),
          plot = gr, width = 10.5, height = 7.0, units = "in"
        )

        plots_list[[ paste0(j, "_vs_", i) ]] = gr

      }
    }
  }
  
  
  return(plots_list)
}
