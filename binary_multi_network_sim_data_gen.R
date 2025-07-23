

#######
# Function: binary_multi_network_sim_data_gen
#
# Description: Generates a 2-layer multilayer network predictor and a binary response 
#              for N observations
#
# Parameters -
# V: Number of nodes for the network on each layer
# N: Number of observations and multilayer networks
# R.star: True latent dimension
# pi.nz: True probability of node activation
# mu.tr: The value of the true location parameter
# seed: The value for the random seed
#
# Returns - 
# y: The binary response (of size N)
# Xmat1: The upper-triangular portion of the network on the first layer for each observation
#        (of dimension N x V(V-1)/2)
# Xmat2: The upper-triangular portion of the network on the second layer for each observation
#        (of dimension N x V(V-1)/2)
#######

binary_multi_network_sim_data_gen <- function(V, N, R.star, pi.nz, mu.tr, seed) {
  
  set.seed(seed)
  
  L <- 2         ## number of network layers
  Q <- V*(V-1)/2 # number of upper-triangular edges on each layer
  
  
  
  ## data generation
  
  Xmat1 <- matrix(rnorm(N*Q),N,Q)  ##upper triangular part of the first undirected network
  Xmat2 <- matrix(rnorm(N*Q),N,Q)  ##upper triangular part of the second undirected network
  
  
  nz.index <- sample(1:V,pi.nz*V,replace=F) ## sampling indices of the influential nodes
  
  # Sampling the latent effects for each layer
  u.tr  <- matrix(0,V,L*R.star)
  
  u1.tr <- array(0,dim=c(V,R.star))
  u2.tr <- array(0,dim=c(V,R.star))
  
  for(v in 1:length(nz.index)){
    
    u.tr[nz.index[v],] <- rnorm(L*R.star)
  }
  
  u1.tr <- u.tr[, c(1:R.star)]
  u2.tr <- u.tr[, c((R.star + 1):(L*R.star))]

  
  
  ### network coefficient
  B1.tr <- tcrossprod(u1.tr) #u1.tr%*%t(u1.tr)
  B2.tr <- tcrossprod(u2.tr) #u2.tr%*%t(u2.tr)
  
  b1.tr <- upperTriangle(B1.tr,diag=FALSE,byrow=T)
  b2.tr <- upperTriangle(B2.tr,diag=FALSE,byrow=T)
  
  
  ## response
  p_y <- exp(mu.tr + Xmat1%*%b1.tr + Xmat2%*%b2.tr) / (1 + exp(mu.tr + Xmat1%*%b1.tr + Xmat2%*%b2.tr)) 
  
  
  y = rep(NA, N)
  
  for (i in 1:N) {
    
    y[i] = rbinom(1, 1, p_y[i])
    
  }
  
  return(list(y = y, Xmat1 = Xmat1, Xmat2 = Xmat2, b1.tr = b1.tr, b2.tr = b2.tr,
              nz.index = nz.index))
  
  
}

# Example function call

#data = binary_multi_network_sim_data_gen(V = 20, N = 250, R.star = 4, pi.nz = 0.2,
#                                          mu.tr = 0.5, seed = 1235)
