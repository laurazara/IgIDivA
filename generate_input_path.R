library(data.table)

generate_input_file = function(input_path){
  files = dir(input_path)
  files=files[which(grepl('.txt',files))]
  ga_files = which(grepl('Grouped Alignment_nt_',files))
  hsim_files = which(grepl('highly_sim_all_clonotypes',files))
  id = strsplit(files,'\\_|\\.')
  id = unlist(lapply(id,function(x){return(x[length(x)-1])}))
  id_ga = id[ga_files]
  id_hsim = id[hsim_files]
  ind_ga = order(id_ga)
  ind_hsim = order(id_hsim)
  id_ga = id_ga[ind_ga]
  ga_files = files[ga_files[ind_ga]]
  id_hsim = id_hsim[ind_hsim]
  hsim_files = files[hsim_files[ind_hsim]]
  if (length(id_ga)==length(id_hsim) & all(id_ga==id_hsim)){
    input_table = data.table(highly_sim_clonos_file=paste0(input_path,'/',hsim_files),
                             grouped_alignment_file=paste0(input_path,'/',ga_files),
                             sample_id=id_hsim,
                             stringsAsFactors=FALSE)
    write.table(input_table,"input_files.txt", sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE,quote = FALSE)
    base::message(paste0('File path: ', getwd(), '/', "input_files.txt"))
    return(input_table)
  }
  else{
    base::message("ids in 'highly_sim_all_clonotypes' but not in 'Grouped Alignment_nt':\n")
    print(setdiff(id_hsim,id_ga))
    base::message("ids in 'Grouped Alignment_nt' but not in 'highly_sim_all_clonotypes':\n")
    print(setdiff(id_ga,id_hsim))
    stop('problem with the ids of the files')
  }
}



generate_output_path = function(save_path){
    
    if(!file.exists(paste0(getwd(), '/', save_path))){
        dir.create(paste0(getwd(), '/', save_path), showWarnings = FALSE)
    }
    
}
