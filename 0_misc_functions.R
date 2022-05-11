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

get_discretization_step <- function(domain_bounds, n_nodes){
  # Either same number of nodes in every direction, or vector
  n_nodes <- c(n_nodes)
  delta_i <- domain_bounds[, 2] - domain_bounds[, 1]
  delta_i <- delta_i / (n_nodes-1)
  return(c(delta_i))
}

index_to_multiindex <- function(index, n_nodes, dim){
  multi_index <- sapply(index, FUN=function(x){
    return(((x-1) %/% (n_nodes^(seq(dim)-1)))%%n_nodes)
  })
  if(dim==1){
    multi_index <- matrix(multi_index, nrow=1)
  }
  return(multi_index)
}

multiindex_to_index <- function(multi_index, n_nodes, dim){
  index <- as.vector(1 + (n_nodes^(seq(dim)-1))%*% multi_index)
  return(index)
}

multiindex_to_coord <- function(multi_index, delta_i, lower_bound){
  coord <- round(delta_i*multi_index+lower_bound, 16)
  return(coord)
}

coord_to_multiindex <- function(coord, delta_i, lower_bound){
  multi_index <- round((coord-lower_bound)/delta_i, 16)
  return(multi_index)
}

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

basis_fun_value <- function(location_multi_index, n_nodes, dim, delta_i){
  n <- ncol(location_multi_index)
  i <- 1
  for(i in seq(n)){
    lower_ind <- floor(location_multi_index[, i])
    upper_ind <- ceiling(location_multi_index[, i])
    if(dim==1){
      multi_index_neighbours_grid <- unname(cbind(lower_ind, upper_ind))
      if(lower_ind==upper_ind){
        multi_index_neighbours_grid<-as.matrix(multi_index_neighbours_grid[, 1])
      }
    }else{
      M <- unname(cbind(lower_ind, upper_ind))
      multi_index_neighbours_grid <- t(
        unique(
          expand.grid(
            lapply(seq(nrow(M)), function(k){
              M[k, ]
            })
          )))
    }
    value <- unname(apply(abs(multi_index_neighbours_grid - 
                                location_multi_index[, i]), 
                          2, max))
    index_neighbours_grid <- multi_index_to_index(multi_index_neighbours_grid, n_nodes, dim)
    if(i==1){
      fun_value <- round(data.frame(index_x = i,
                                    index_f=index_neighbours_grid,
                                    value=1-value), 16)
      fun_value <- unique(fun_value)
    }else{
      Fun_temp <- round(data.frame(index_x = i,
                                   index_f=index_neighbours_grid,
                                   value=1-value), 16)
      Fun_temp <- unique(Fun_temp)
      fun_value <- rbind(fun_value, Fun_temp)
    }
  }
  fun_value <- round(fun_value, 14)
  fun_value <- fun_value[fun_value$value!=0, ]
  return(fun_value)
}

eval_gp <- function(fun_value, weight){
  res <- fun_value %>% 
    group_by(index_x)%>%
    summarise(value=sum(value*weight[index_f]))%>%
    dplyr::select(value)
  return(res$value)
}

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

