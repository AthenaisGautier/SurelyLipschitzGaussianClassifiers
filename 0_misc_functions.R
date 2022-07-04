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


#### predict_probability_fullcalc
# A function that takes a (or several) new directions, values for epsilon (weights), the domain bounds and number of nodes, and predicts the probability at considered new_loc by computing everything
####
predict_probability_fullcalc <- function(new_loc, weights, dim,
                                         domain_bounds,
                                         n_nodes){
  delta_i <- get_delta_i(domain_bounds, n_nodes)
  location_multi_index <- coord_to_multi_index(new_loc, delta_i, domain_bounds[, 1])
  basis_fun_value <- basis_fun_value(location_multi_index=location_multi_index, 
                                     n_nodes, dim, delta_i)
  
  GP_values <- sapply(seq(nrow(weights)), function(i){
    eval_gp(basis_fun_value, weights[i, ])
  })
  sigmoid_values <- 1/(1+exp(-GP_values)) 
  df <- data.frame(x=t(new_loc), GP=GP_values, sigmoid=sigmoid_values)
  return(df)
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



