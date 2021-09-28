generate_input_file = function(input_path,input_file_name=""){
  files = dir(input_path)
  if (input_file_name==""){
    input_file_name = paste0(dirname(input_path),'/input_files.txt')
  }
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
    input_table = data.frame(highly_sim_clonos_file=paste0(input_path,'/',hsim_files),
                             grouped_alignment_file=paste0(input_path,'/',ga_files),
                             sample_id=id_hsim,
                             stringsAsFactors=FALSE)
    write.table(input_table,input_file_name, sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE,quote = FALSE)
    base::message(paste0('File path: ',input_file_name))
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
