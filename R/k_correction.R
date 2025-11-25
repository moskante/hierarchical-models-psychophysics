compute_k = function(p, gamma, lambda){
  
  k <- qnorm( (p - gamma)/(1 - gamma - lambda) )
  
  return(k)
}

#example
compute_k(0.5, 0.027, 0.024)
