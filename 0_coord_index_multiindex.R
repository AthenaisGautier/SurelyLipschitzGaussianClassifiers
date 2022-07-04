#### index_to_multiindex
# A function that maps indices (1D) to multi-indices (dim-D)
####
index_to_multiindex <- function(index, n_nodes, dim){
  multi_index <- sapply(index, FUN=function(x){
    return(((x-1) %/% (n_nodes^(seq(dim)-1)))%%n_nodes)
  })
  if(dim==1){
    multi_index <- matrix(multi_index, nrow=1)
  }
  return(multi_index)
}

#### multiindex_to_index
# A function that maps multi-indices (dim-D) to indices (1D)
####
multiindex_to_index <- function(multi_index, n_nodes, dim){
  index <- as.vector(1 + (n_nodes^(seq(dim)-1))%*% multi_index)
  return(index)
}

#### multiindex_to_coord
# A function that maps multi-indices (dim-D) to coordinates (dim-D)
####
multiindex_to_coord <- function(multi_index, delta_i, lower_bound){
  coord <- round(delta_i*multi_index+lower_bound, 16)
  return(coord)
}

#### coord_to_multiindex
# A function that maps coordinates (dim-D) to multi-indices (dim-D)
####
coord_to_multiindex <- function(coord, delta_i, lower_bound){
  multi_index <- round((coord-lower_bound)/delta_i, 16)
  return(multi_index)
}
