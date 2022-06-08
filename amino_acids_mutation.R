library(RGenetics)
library(data.table)

get_nodes_aa_muts = function(nodes,germline){
#the germline has to start at position 1
germline = as.character(germline)
nt_muts = unique(do.call(c,strsplit(nodes$Mutations[nodes$suppl_mu>=0],' ')))
n_muts = length(nt_muts)
n_nodes = nrow(nodes)
mut_table = data.table(nt_muts=nt_muts,from=character(n_muts),to=character(n_muts),pos=integer(n_muts),aa_pos=integer(n_muts),codon_start=integer(n_muts),aa_origin=character(n_muts),stringsAsFactors = FALSE)

for (i in 1:n_muts){
  mut = strsplit(nt_muts[i],"")[[1]]
  l = length(mut)
  pos = as.integer(paste(mut[2:(l-1)],collapse=""))
  aa_pos = floor((pos-1)/3+1)
  codon_ind = 3*(aa_pos-1)+1 
  codon_o = germline[codon_ind:(codon_ind+2)]
  codon_o = paste(codon_o,collapse="")
  aa_o = codonToAAone(codon_o)
  
  mut_table$from[i] = mut[1]
  mut_table$to[i] = mut[l]
  mut_table$pos[i] = pos
  mut_table$aa_pos[i] = aa_pos 
  mut_table$codon_start[i] = codon_ind
  mut_table$aa_origin[i] = aa_o
}
ind_pos = which(nodes$suppl_mu>0)
ind_pos = append(which(nodes$suppl_mu==0),ind_pos)
nodes$Mutations_aa = character(n_nodes)
nodes$count_mut_aa = numeric(n_nodes)
germline_main_variant = germline
for (i in ind_pos){
  germline_mut = germline_main_variant
  muts_curr = strsplit(nodes$Mutations[i]," ")[[1]]
  for (m in muts_curr){
    ind = which(mut_table$nt_muts==m)
    germline_mut[mut_table$pos[ind]] = mut_table$to[ind]
  }
  nb_mut_curr = length(muts_curr)
  muts_aa_curr = matrix("",nrow=1,ncol=nb_mut_curr)
  if (nb_mut_curr>0){
    for (j in 1:nb_mut_curr){
      m = muts_curr[j]
      ind = which(mut_table$nt_muts==m)
      codon_e = germline_mut[mut_table$codon_start[ind]:(mut_table$codon_start[ind]+2)]
      if (length(setdiff(codon_e,c('a','c','g','t')))==0){
        codon_e = paste(codon_e,collapse="")
        aa_e = codonToAAone(codon_e)
        if (aa_e!=mut_table$aa_origin[ind]){
          aa_mut= sprintf("%s%i%s\n",mut_table$aa_origin[ind],mut_table$aa_pos[ind],aa_e)
          muts_aa_curr[muts_aa_curr==aa_mut] = ""
          muts_aa_curr[j] = aa_mut
        }
      }
    }
  }
  col2collapse = muts_aa_curr!=""
  muts_aa_curr = muts_aa_curr[col2collapse]
  if (nodes$suppl_mu[i]==0){
    germline_main_variant = germline_mut
    muts_aa_main_variant = muts_aa_curr
  }
  else {
    muts_aa_curr = setdiff(muts_aa_curr,muts_aa_main_variant)
  }
  nodes$count_mut_aa[i] = sum(col2collapse)
  nodes$Mutations_aa[i] = paste(muts_aa_curr,collapse="")
}
nodes$count_mut_aa[nodes$suppl_mu<0] = -1
return(nodes)
}

get_aa_muts_weight = function(nodes,with_main_nt_variant=FALSE){
  if (with_main_nt_variant){
    aa_muts = unique(do.call(c,strsplit(nodes$Mutations_aa[nodes$suppl_mu==0],'\n')))
  }
  else {
    aa_muts = unique(do.call(c,strsplit(nodes$Mutations_aa[nodes$suppl_mu>0],'\n')))
  }
  n_aa_muts = length(aa_muts)
  if (n_aa_muts>0){
    aa_muts_weight = data.table(aa_muts=aa_muts, seqs=integer(n_aa_muts))
    for (i in 1:n_aa_muts){
      aa_muts_weight$seqs[i] = sum(nodes$N[grepl(aa_muts[i],nodes$Mutations_aa,fixed=TRUE)])
    }
    aa_muts_weight = aa_muts_weight[
      with(aa_muts_weight, order(seqs, decreasing = TRUE)),
      ]
  }
  else {aa_muts_weight = data.table(aa_muts=character(0), seqs=integer(0))}
  return(aa_muts_weight)
}
