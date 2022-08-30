#### prepare_data
# A function that takes a messy dataset (with duplicates in the x) and returns a 
# Clean table of location, successes and tries. Uses dplyr 
####
prepare_data <- function(data, name_index, req.dplyr=TRUE){
  name_index <- c(name_index)
  if(!all(data$class %in% c(-1, 0, 1))){
    warning("The `class` variable should take value 1 (for success) or 0/-1.")
  }
  if(req.dplyr){
    require(dplyr)
    data_clean <- data%>%
      group_by_at(name_index)%>%
      summarise(multiplicity=n(),
                success=sum(class==1),
                .groups="keep")%>%
      data.frame()
  }else{
    error("Was not implemented yet without dplyr")
  }
  return(data_clean)
}

#### get_discretization_step
# A function that takes the domain bound in each direction, and gives the step size between nodes along each direction
####
get_discretization_step <- function(domain_bounds, n_nodes){
  # Either same number of nodes in every direction, or vector
  n_nodes <- c(n_nodes)
  delta_i <- domain_bounds[, 2] - domain_bounds[, 1]
  delta_i <- delta_i / (n_nodes-1)
  return(c(delta_i))
}

#### basis_fun_value
# A function that takes a location_multi_index, the indices of the nodes,
# and other data, and returns a matrix of the basis functions centered at the nodes evaluated at the new locations
####
basis_fun_value <- function(location_multi_index, multi_index_nodes, n_nodes, dim, delta_i){
  n <- nrow(location_multi_index)
  return(t(sapply(seq(n), function(k){
    x <- location_multi_index[k, ]
    round(pmax(0, 1-pmax(t(abs(multi_index_nodes-x)))), 15)
  })))
}


#### field_predict
# A function that takes a (or several) new locations, values for epsilon (weights), the domain bounds and number of nodes, and predicts the probability at considered new_loc by computing everything
####

field_predict <- function(new_coord, weights, name_index, domain_bounds, dim,
                          n_nodes){
  if(is.null(dim(weights))){
    weights <- t(as.matrix(weights))
  }
  # Useful quantities regarding the nodes we consider
  delta_i <- get_discretization_step(domain_bounds, n_nodes)
  multi_index_nodes <- index_to_multiindex(seq(1, n_nodes^dim), n_nodes, dim)
  coord_nodes <- t(multiindex_to_coord(multi_index_nodes, delta_i, domain_bounds[, 1]))
  colnames(coord_nodes) <- name_index
  new_coord_loc <- new_coord[, name_index, drop=FALSE]
  location_multi_index <- coord_to_multiindex(coord=new_coord_loc, 
                                              delta_i, lower_bound = domain_bounds[, 1])
  fun_value <- basis_fun_value(location_multi_index=location_multi_index,
                               multi_index_nodes=multi_index_nodes,
                               n_nodes, dim, delta_i)
  GP_val <- fun_value%*%  t(weights)
  res <- data.frame(ix = seq(nrow(new_coord)), 
                    GP=c(GP_val), 
                    real=c(sapply(seq(nrow(weights)), FUN=function(k){
                      rep(k, nrow(new_coord))})))
  res <- cbind(new_coord[res$ix, , drop=FALSE], res[, -c(1)])
  res$prob <- 1/(1+exp(-res$GP))
  return(res)
}

#### constraints_matrix
# A function that creates the constraint matrix Fmat and vector g_i such that Fmat %*% weight <= g_i
####
constraints_matrix <- function(multi_index_nodes, dim, 
                               n_nodes, 
                               C_lip, delta_i){
  #We need to identify which nodes are neihbours to one another
  following_mat <- diag(dim)
  for(i in seq(n_nodes^dim)){
    current_node <- multi_index_nodes[, i]
    followers <- current_node+following_mat
    admissible <- colSums(followers < n_nodes)==dim
    n_adm <- sum(admissible)
    ind_admissible <- multiindex_to_index(followers[, admissible], n_nodes, dim)
    
    if(i==1){
      F_mat <- matrix(0, nrow=n_adm, ncol=n_nodes^dim)
      F_mat[, i] <- -1
      F_mat[seq(n_adm), ind_admissible] <- diag(n_adm)
      g_i <- C_lip * delta_i[admissible]
    }else{
      F_temp <- matrix(0, nrow=n_adm, ncol=n_nodes^dim)
      F_temp[, i] <- -1
      F_temp[seq(n_adm), ind_admissible] <- diag(n_adm)
      F_mat <- rbind(F_mat, F_temp)
      g_temp <- C_lip * delta_i[admissible]
      g_i <- c(g_i, g_temp)
    }
  }
  return(list(Fmat=rbind(F_mat, -F_mat), g_i=c(g_i, g_i)))
}



