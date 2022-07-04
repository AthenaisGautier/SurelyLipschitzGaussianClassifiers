ibrary(numDeriv)


#### prior_lengthscale
# A function that gives the prior over the lengthscale (or its log) and its derivative 
# (so far we use an log-normal 0,3 and we use numerical derivatives)
####
prior_lengthscale <- function(theta, isotropy=TRUE, log=TRUE, dim=1){
  if(!isotropy){
    if(length(theta)!=dim){
      stop("The vector `theta` has the wrong size (should match `dim`)")
    }
  }  
  require(numDeriv)
  res <- dlnorm(theta, meanlog = 0, sdlog = 3, log = log)
  if(log){
    res <- sum(res)
  }else{
    res <- prod(res)
  }
  derive <- grad(func=function(x){res <- dlnorm(x, meanlog = 0, sdlog = 3, log = log)
    if(log){
      res <- sum(res)
    }else{
      res <- prod(res)
    }
  return(res)})
  return(list(value=res, grad=derive))
}

