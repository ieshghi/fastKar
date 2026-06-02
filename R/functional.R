#use ABC model to score enhancer interactions with promoters. 
#returns matrix of dimension N_enhancer x N_promoter
#Each promoter column sums to 1, so the entries of the matrix are probabilities of an enhancer activating a given promoter
abc_from_tiles = function(enhancers, promoters, contact.gm ,enhancer_score_col = "score",promoter_name_col = "gene",eps = 1e-12) {
  tiles = contact.gm$gr
  contacts = contact.gm$dat
  n_tiles = length(tiles)
  n_enh = length(enhancers)
  n_prom = length(promoters)
  C = sparseMatrix(
    i = c(contacts$i,contacts$j),
    j = c(contacts$j,contacts$i),
    x = rep(contacts$value,2),
    dims = c(n_tiles, n_tiles)
  )
  # tile -> enhancer incidence
  te = gr.findoverlaps(tiles, enhancers)
  T_E = sparseMatrix(i = te$query.id,j = te$subject.id,x = 1,dims = c(n_tiles, n_enh))
  # tile -> promoter incidence
  tp = gr.findoverlaps(tiles, promoters)
  T_P = sparseMatrix(i = tp$query.id,j = tp$subject.id,x = 1,dims = c(n_tiles, n_prom))
  # enhancer x promoter contact mass
  EP_contact = t(T_E) %*% C %*% T_P
  # weight each enhancer by activity score
  enh_score = mcols(enhancers)[[enhancer_score_col]]
  EP_raw = EP_contact * enh_score
  # normalize within each promoter
  denom = colSums(EP_raw) + eps
  EP_abc = sweep(as.matrix(EP_raw), 2, denom, "/")
  rownames(EP_abc) = seq_len(n_enh)
  colnames(EP_abc) = mcols(promoters)[[promoter_name_col]]
  return(EP_abc)
}
