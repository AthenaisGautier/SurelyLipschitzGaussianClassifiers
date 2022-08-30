library(rlang)

#### prepare_optim
# A function that prepares the data and code
# To perform an optim (finding hyperparameters)
####
prepare_optim <- function(data, name_index, domain_bounds, dim,
                          n_nodes, C_lip, kernel_choice, param_lengthscale,
                          mean_val){
  
  data_clean <- prepare_data(data, name_index)
  require(kergp)
  if(kernel_choice=="Mat52"){
    cov <- kMatern(d=dim,
                   nu="5/2")
  }else{
    if(kernel_choice=="Mat32"){
      cov <- kMatern(d=dim,
                     nu="3/2")
    }else{
      if(kernel_choice=="Exp"){
        cov <- kExp(d=dim)
      }else{
        if(kernel_choice=="Gauss"){
          cov <- kGauss(d=dim)
        }
      }
    }
  }
  inputNames(cov) <- name_index
  coef(cov) <- c(param_lengthscale, 1)
  range_domain <- apply(domain_bounds, 1, diff)
  coefUpper(cov) <- c(range_domain, 1)
  coefLower(cov) <- c(range_domain/100, 1)
  
  # Useful quantities regarding the nodes we consider
  delta_i <- get_discretization_step(domain_bounds, n_nodes)
  multi_index_nodes <- index_to_multiindex(seq(1, n_nodes^dim), n_nodes, dim)
  coord_nodes <- t(multiindex_to_coord(multi_index_nodes, delta_i, domain_bounds[, 1]))
  colnames(coord_nodes) <- name_index
  
  #Get constraint matrix
  cstMat <- constraints_matrix(multi_index_nodes, dim, 
                               n_nodes, 
                               C_lip, delta_i)
  
  # # Get covariance of nodes
  # C <- covMat(cov, coord_nodes)
  #Check the mean values provided at nodes
  if(is.function(mean_val)){
    mean_val <- apply(coord_nodes, 1, mean_val)
  }
  mean_val <- c(mean_val)
  if(length(mean_val)==1){
    mean_val <- rep(mean_val, n_nodes^dim)
  }
  if(length(mean_val)!=n_nodes^dim){
    stop(paste0("You provided a mean vector's size does not fit. (Length is ", 
                length(mean_val), " expected ", n_nodes^dim))
  }
  
  #Get functions value at samples 
  location_data <- data_clean[, name_index, drop=FALSE]
  location_multi_index <- coord_to_multiindex(coord=location_data, 
                                              delta_i, lower_bound = domain_bounds[, 1])
  fun_value <- basis_fun_value(location_multi_index=location_multi_index,
                               multi_index_nodes=multi_index_nodes,
                               n_nodes, dim, delta_i)
  
  
  return(list(fun_value=fun_value,
              data_clean=data_clean,
              coord_nodes=coord_nodes,
              cov=cov,
              mean_prior=mean_val,
              cstMat=cstMat))
}

#### EvaluatePostAndItsderiv
# A function that  evaluates the negative log posterior (up to a constant) and its
# derivatives. Useful when one need to optimise both the lengthscales and weights
####
EvaluatePostAndItsderivForBoth <- function(weights, dim,
                                           lengthscale=NULL, kernel,
                                           mean_prior, coord_nodes,
                                           successes, tries,
                                           fun_value,
                                           comp_lengthscale=TRUE,
                                           sd_logtheta = 3,
                                           comp_grad_weights=TRUE,
                                           comp_grad_len=TRUE,
                                           Sigma=NULL, 
                                           invSigma=NULL,
                                           logdetSigma=NULL){
  centered_weights <- weights - mean_prior
  
  if(is.null(Sigma)){
    # We did not provide the covariance matrix Sigma of the weights
    coef(kernel) <-c(lengthscale, 1)
    Sigma <- covMat(kernel, coord_nodes)
  }
  if(is.null(invSigma) | is.null(logdetSigma)){
    # We did not provide the inverse or log-determinant of Sigma
    cholSigma <- chol(Sigma/2+t(Sigma)/2+1e-9*diag(nrow(Sigma)))
    if(is.null(invSigma)){
      invSigma <- chol2inv(cholSigma)
    }
    if(is.null(logdetSigma)){
      logdetSigma <- sum(log(diag(cholSigma)))*2
    }
  }
  p <- nrow(coord_nodes)
  prior_term_1 <- invSigma %*% centered_weights
  prior_term_2 <- t(centered_weights) %*% prior_term_1
  GP_value <- as.vector(fun_value %*% weights)
  tempexp <- exp(-GP_value)
  sigmoid <- 1/(1+tempexp)
  
  res <-  logdetSigma/2 + prior_term_2/ 2 - # Term from the prior on weights
    sum(successes*log(sigmoid)) - sum((tries-successes)*log(1-sigmoid)) # Term from the likelihood
  if(comp_lengthscale){
    # We compute the likelihood WITH contribution from the lengthscale
    if(!is.null(lengthscale)){
      res <- res + sum(log(lengthscale)) + sum(log(lengthscale)^2)/(2*sd_logtheta^2)
    }else{
      stop("Must provide the value `lengthscale` if `comp_lengthscale` is TRUE")
    }
  }
  if(comp_grad_weights){
    grad_weights <- prior_term_1 -
      colSums(successes*(1-sigmoid)*fun_value) +
      colSums((tries-successes)*sigmoid*fun_value)
  }else{
    grad_weights <- NULL
  }
  if(comp_grad_len){
    grad_theta  <- sapply(seq(dim), FUN=function(k){
      sum(diag(invSigma%*% attr(Sigma, "grad")[,,k])) /2 }) +
      sapply(seq(dim), FUN=function(k){
        t(centered_weights) %*%
          (- invSigma %*% attr(Sigma, "grad")[,,k]%*% invSigma) %*% 
          centered_weights/2  })+
      ((sd_logtheta^2+log(lengthscale))/(sd_logtheta^2*lengthscale))
  }else{
    grad_theta  <- NULL
  }
  return(list(value=res, 
              grad_weights=as.vector(grad_weights), 
              grad_theta=as.vector(grad_theta)))
}


#### run_optim_twosteps
# A function that performs an optim with Gibbs steps
# Currently does not work...
####

run_optim_twosteps <- function(prepared_optim, sd_logtheta=3,
                               weight_start, lengthscale_start, 
                               xrel, frel, maxstep, print_level=0){
  current_weights <- weight_start
  current_lengthscale <- lengthscale_start
  env_optim <- new.env(parent = empty_env())
  cleaned_data <- prepared_optim$data_clean
  cov <- prepared_optim$cov
  step <- 1
  converged <- FALSE
  track<-c()
  for(step in seq(maxstep)){
    if(converged){
      break
    }
    if(print_level>=1){
      print(paste0("Step ", step))
    }
    # Step on the weights
    coef(cov) <- c(current_lengthscale, 1)
    Sigma <- covMat(cov, 
                    prepared_optim$coord_nodes)
    invSigma <- solve(Sigma)
    logdetSigma<-log(det(Sigma))
    res_opt <- constrOptim(theta=current_weights, 
                           f=function(x){
                             res <-EvaluatePostAndItsderivForBoth(weights=x, 
                                                                  dim=dim,
                                                                  lengthscale=current_lengthscale, 
                                                                  kernel=cov,
                                                                  mean_prior=prepared_optim$mean_prior,
                                                                  coord_nodes=prepared_optim$coord_nodes,
                                                                  successes=cleaned_data$success, 
                                                                  tries=cleaned_data$multiplicity,
                                                                  fun_value=prepared_optim$fun_value,
                                                                  comp_lengthscale =TRUE,
                                                                  sd_logtheta = sd_logtheta,
                                                                  comp_grad_weights=TRUE,
                                                                  comp_grad_len=FALSE,
                                                                  Sigma=Sigma, 
                                                                  invSigma=invSigma,
                                                                  logdetSigma=logdetSigma) 
                             
                             env_optim$fn <- res$value
                             env_optim$grad <- res$grad_weights
                             return(res$value)
                           },
                           grad=function(x){
                             return(env_optim$grad)
                           },
                           ui=-prepared_optim$cstMat$Fmat,
                           ci=-prepared_optim$cstMat$g_i)
    new_weights <- res_opt$par
    if(step==1){
      old_p_weights <- Inf
      old_p_lengthscale <- Inf
    }
    track <- c(track, res_opt$value[1])
    current_p_weights <- res_opt$value[1]
    increment_weight <- max(abs((new_weights-current_weights)/current_weights))
    increment_p_weights <- abs((old_p_weights-current_p_weights)/old_p_weights)
    old_p_weights <- current_p_weights
    current_weights <- new_weights
    #Step on the hyperparameters
    require(rgenoud)
    res_opt <- genoud(fn=function(x){
      res <- EvaluatePostAndItsderivForBoth(weights=current_weights, dim=dim,
                                            lengthscale=x, kernel=cov,
                                            mean_prior=prepared_optim$mean_prior, 
                                            coord_nodes=prepared_optim$coord_nodes,
                                            successes=cleaned_data$success, 
                                            tries=cleaned_data$multiplicity,
                                            fun_value=prepared_optim$fun_value,
                                            Sigma=NULL,invSigma=NULL, 
                                            comp_lengthscale =TRUE,
                                            sd_logtheta = sd_logtheta,
                                            comp_grad_weights=FALSE,
                                            comp_grad_len=TRUE)
      env_optim$fn <- res$value
      env_optim$grad <- res$grad_theta
      return(res$value)
    },
    gr=function(x){
      return(env_optim$grad)
    },
    max=FALSE, gradient.check = FALSE, nvars=dim,
    Domains=matrix(c(coefLower(cov)[1:dim], coefUpper(cov)[1:dim]), ncol=2),
    BFGS=TRUE, boundary.enforcement = 2, print.level = 0)
    new_lengthscale <- res_opt$par
    increment_lengthscale <- max(abs((new_lengthscale-current_lengthscale)/current_lengthscale))
    track <- c(track, res_opt$value[1])
    current_p_lengthscale <- res_opt$value[1]
    increment_p_lengthscale <- abs((old_p_lengthscale-current_p_lengthscale)/old_p_lengthscale)
    old_p_lengthscale <- current_p_lengthscale
    current_lengthscale <- new_lengthscale
    if(step >1){
      if(increment_p_lengthscale <frel & increment_p_weights < frel){
        print("Changes in the objective function are below the relative threshold")
        converged<- TRUE
      }
      if(increment_lengthscale <xrel & increment_weight < xrel){
        print("Changes in the parameters are below the relative threshold")
        converged <- TRUE
      }
    }
    
  }
  if(!converged){
    warning("The algorithm stopped because the maximum number is reached")
  }
  track <- data.frame(type=c("weight", "lengthscale"), value=track)
  track$step <- (round((seq(nrow(track))-1)/2))+1
  track$i <- seq(nrow(track))
  return(list(lengthscale=current_lengthscale, weights=current_weights, track=track))
}


#### run_optim_fixed_len
# A function that performs an optim with a given lengthscale
####

run_optim_fixed_len <- function(prepared_optim, sd_logtheta=3,
                                lengthscale_candidates){
  env_optim <- new.env(parent = empty_env())
  cleaned_data <- prepared_optim$data_clean
  cov <- prepared_optim$cov
  res_opt <- constrOptim(theta=weight_start, 
                         f=function(x){
                           res <-EvaluatePostAndItsderivForBoth(weights=x, 
                                                                dim=dim,
                                                                lengthscale=lengthscale_candidates, 
                                                                kernel=cov,
                                                                mean_prior=prepared_optim$mean_prior,
                                                                coord_nodes=prepared_optim$coord_nodes,
                                                                successes=cleaned_data$success, 
                                                                tries=cleaned_data$multiplicity,
                                                                fun_value=prepared_optim$fun_value,
                                                                comp_lengthscale =TRUE,
                                                                sd_logtheta = sd_logtheta,
                                                                comp_grad_weights=TRUE,
                                                                comp_grad_len=FALSE,
                                                                Sigma=Sigma, 
                                                                invSigma=invSigma,
                                                                logdetSigma=logdetSigma) 
                           
                           env_optim$fn <- res$value
                           env_optim$grad <- res$grad_weights
                           return(res$value)
                         },
                         grad=function(x){
                           return(env_optim$grad)
                         },
                         ui=-prepared_optim$cstMat$Fmat,
                         ci=-prepared_optim$cstMat$g_i)
  
  return(list(lengthscale=lengthscale_candidates, weights=res_opt$par, value=res_opt$value))
}
