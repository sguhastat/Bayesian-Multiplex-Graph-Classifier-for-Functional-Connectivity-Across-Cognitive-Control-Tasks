library(pscl)
library(mvtnorm)
library(gdata)
library(LaplacesDemon)
library(MASS)
library(GIGrvg)
library(mvtnorm)
library(MCMCpack)
library(BayesLogit)

library(pROC)


binary_multilayer_network_model <- function (y, Xmat1, Xmat2, W = matrix(0, nrow = length(y), ncol = 2), 
                                                            b1.tr = NULL, b2.tr = NULL,
                                                            niter = 10000, 
                                                            n_burnin = 0,
                                                            nthinning = 1,
                                                            seed = 1234, test_size = 0.2,
                                                            L = 2, R = 6, a = 1, b = 1, pi = a / (a+b), 
                                                            M_scale = 1, normalize_M = FALSE) {
  
  ## model fitting
  set.seed(seed)
  
  # Number of upper triangular entries in each network
  Q = ncol(Xmat1)
  
  L = 2
  p = ncol(W)
  
  # Getting number of vertices in each network based on the number of upper triangular entries
  V = (1 + sqrt(1 + 8 * Q)) / 2
  
  
  # Splitting data into train and test set
  test_idxs = sample(c(1:length(y)), size = test_size * length(y), replace = FALSE)
  
  Xmat1_test = Xmat1[test_idxs, ]
  Xmat2_test = Xmat2[test_idxs, ]
  W_test = W[test_idxs, ]
  
  y_test = y[test_idxs]
  
  
  Xmat1 = Xmat1[-test_idxs, ]
  Xmat2 = Xmat2[-test_idxs, ]
  W = W[-test_idxs, ]
  
  y = y[-test_idxs]
  
  N = length(y)
  
  
  
  a.tau <- 1       
  b.tau <- 1
  r <- 1
  delta <- 1 
  #a <- 1
  #b <- 1
  a.wish <- L * R + 1
  S.scale <- (1 / (M_scale * a.wish)) * diag(rep(1,L * R))
  eta = 1.01
  
  #### Hyperparameters for lambda probability
  a.kap <- c()
  b.kap <- c()
  c.kap <- c()
  
  for(r in 1:R){
    a.kap[r] <- 1
    b.kap[r] <- r^(eta)
    c.kap[r] <- 1
  }
  
  #grid.index <- expand.grid(1:V,1:V)
  ####################################################################################
  #### storing MCMC iterates of the parameters
  ####################################################################################
  

  mu.store  <- numeric()
  
  gamma.store <- matrix(NA, niter, p)
  
  pi.store <- numeric()

  
  xi.store <- list() ## Same as matrix M

  
  rn.gen.store <- matrix(NA,niter,V)
  

  
  kappa1.store  <- matrix(NA,niter,R)
  kappa2.store  <- matrix(NA,niter,R)
  

  pi1.kap.store <- array(NA,dim=c(niter,R,3))
  pi2.kap.store <- array(NA,dim=c(niter,R,3))
  
  ulist.store <- list()  ## each iterate will store a V\times L*R matrix
  ulist1.store <- list() ## each iterate will store a V\times R matrix
  ulist2.store <- list() ## each iterate will store a V\times R matrix
  
  
  # Add storage for omega, matrix of size niter by N
  omega.store <- matrix(NA, niter, N)
  
  #### Gibbs Sampler
  #### Initial Values

  mu <- 0.5

  pi1.kap <- matrix(NA,R,3)
  pi2.kap <- matrix(NA,R,3)
  

  kappa1 <- numeric()
  kappa2 <- numeric()
  
  for(r in 1:R){

    pi1.kap[r,] <- c(rdirichlet(1,c(a.kap[r],b.kap[r],c.kap[r])))
    pi2.kap[r,] <- c(rdirichlet(1,c(a.kap[r],b.kap[r],c.kap[r])))
    
    #kappa[r]  <- sample(c(1,0,-1),1,pi.kap[r,],replace=T)
    kappa1[r] <- sample(c(1,0,-1),1,pi1.kap[r,],replace=T)
    kappa2[r] <- sample(c(1,0,-1),1,pi2.kap[r,],replace=T)
    
  }
  

  
  M <- rinvwishart(a.wish,S.scale)

  #### Each of the 'niter' no. of matrices in the list represents u1,u2,u3,...,uV.(Each of the u's is of dimension R.)
  
  # Initialize omega
  omega = rpg(num = N, h = 1, z = 0.0)
  
  
  ulist   <- matrix(rnorm(V*L*R,0,1),nrow = V,ncol= L*R)
  ulist1  <- ulist[, c(1:R)]
  ulist2  <- ulist[, c((R+1):(L*R))]
  
  gamma = c(rnorm(p))
  

  Byatha1 <- diag(kappa1)
  Byatha2 <- diag(kappa2)
  
  

  
  
  # Network matrices for each sample
  new.Xmat1 <- list()
  new.Xmat2 <- list()
  

  
  # (?) Constructing Network Matrices for each sample 
  # (?) Can be done outside loop
  for(nn in 1:N){
    
    # Filling Lower Triangular Portion
    new.X1 <- matrix(0,V,V)
    new.X2 <- matrix(0,V,V)
    
    upperTriangle(new.X1,diag=FALSE,byrow=T) <- Xmat1[nn,]
    upperTriangle(new.X2,diag=FALSE,byrow=T) <- Xmat2[nn,]
    
    new.Xmat1.temp <- t(new.X1)
    new.Xmat2.temp <- t(new.X2)
    
    
    # Filling Upper Triangular Portion
    upperTriangle(new.Xmat1.temp,diag=FALSE,byrow=T) <- Xmat1[nn,]
    upperTriangle(new.Xmat2.temp,diag=FALSE,byrow=T) <- Xmat2[nn,]
    
    
    # Allocating network matrices for nn-th sample
    new.Xmat1[[nn]] <- new.Xmat1.temp
    new.Xmat2[[nn]] <- new.Xmat2.temp
    
  }
  
  
  
  for (i in 1:niter){
    

    
    big_omega = diag(1 / omega)
    big_omega_inv = diag(omega)
    
    k = (y - 0.5) / omega
    
    #### 1.Update tau^2

    coef1 <- upperTriangle(ulist1%*%(t(ulist1) * kappa1),diag=FALSE,byrow=T)
    coef2 <- upperTriangle(ulist2%*%(t(ulist2) * kappa2),diag=FALSE,byrow=T)
    

    mean.fit <- W%*%gamma + Xmat1%*%coef1 + Xmat2%*%coef2
    
    epsilon <- (k - mu - mean.fit)

    

    
    if (i %% 20 == 0 | i == 1) {
      
      print(i)
      
    }
    
    
    #### 2.Update mu

    mu <- rnorm(1, sum(omega * (k - mean.fit)) / (sum(omega)),
                1 / (sum(omega)))
    
    
    #### 2.5 Update gamma
    sigma_gamma <- chol2inv(chol(t(W)%*%big_omega_inv%*%W
                                 + diag(rep(1, p))))
    mu_gamma <- sigma_gamma%*%crossprod(W, c(omega * (k - mu - Xmat1%*%coef1
                                                      - Xmat2%*%coef2)))
    
    gamma <- rmvnorm(1, mean = mu_gamma, sigma = sigma_gamma)
    gamma <- c(gamma)
    
    
    #### 3. Update ulist
    

    
    # Perform joint update of ulist1 and ulist2
    rn.gen <- numeric()
    
    for (v in 1:V) {
     
      ind.v <- c(1:V)[-v]
      
      
      mat1 <- ulist1%*%(t(ulist1) * kappa1)
      X.tilde1 <- matrix(unlist(lapply(1:N,function(nn){new.Xmat1[[nn]][-v,v]})),N,V-1,byrow=T)  
      U.tilde1 <- ulist1[-v,]
      F.tilde1 <- X.tilde1%*%U.tilde1%*%diag(kappa1)
    
      mat2 <- ulist2%*%(t(ulist2) * kappa2)
      X.tilde2 <- matrix(unlist(lapply(1:N,function(nn){new.Xmat2[[nn]][-v,v]})),N,V-1,byrow=T)  
      U.tilde2 <- ulist2[-v,]
      F.tilde2 <- X.tilde2%*%U.tilde2%*%diag(kappa2)
      
      F.tilde = cbind(F.tilde1, F.tilde2)
      
      k.tilde <- (c(k) - 
                  sapply(1:N,function(nn){sum(new.Xmat1[[nn]][ind.v,ind.v]*mat1[ind.v,ind.v])/2}) -  
                  sapply(1:N,function(nn){sum(new.Xmat2[[nn]][ind.v,ind.v]*mat2[ind.v,ind.v])/2}) -
                  W%*%gamma - mu)
      k.tilde = c(k.tilde)
      
      
      log.v1 <- dmvnorm(k.tilde, rep(0,N) ,
                        F.tilde%*%M%*%t(F.tilde) + big_omega,log=T)
      log.v2 <- dmvnorm(k.tilde, rep(0,N), big_omega,log=T)
      prob.v <- pi/((1-pi)*exp(log.v2-log.v1)+pi)
      rn.gen[v] <- rbinom(1,1,prob.v)
  
      if(rn.gen[v]==1){
        Sigma.v <- chol2inv(chol(t(F.tilde)%*%big_omega_inv%*%F.tilde
                                 + chol2inv(chol(M))))
        
        upperTriangle(Sigma.v, diag = FALSE, byrow = TRUE) <- (rnorm(1, mean = 0, sd = 10^(-5)) + 
          upperTriangle(Sigma.v, diag = FALSE, byrow = TRUE))
        lowerTriangle(Sigma.v, diag = FALSE, byrow = FALSE) <- upperTriangle(Sigma.v, 
                                                              diag = FALSE, byrow = TRUE)
  
        mu.v  <- Sigma.v%*%crossprod(F.tilde, c(omega * k.tilde))
  
        ulist[v,] <- rmvnorm(1,mu.v,Sigma.v)
      }else{
        ulist[v,] <- rep(0,L*R)
      }
      
      ulist1[v, ]  <- ulist[v, c(1:R)]
      ulist2[v, ]  <- ulist[v, c((R+1):(L*R))]
    }
    
    coef1 <- upperTriangle(ulist1%*%(t(ulist1) * kappa1),diag=FALSE,byrow=T)
    coef2 <- upperTriangle(ulist2%*%(t(ulist2) * kappa2),diag=FALSE,byrow=T)
    
    #### 8. Update pi (corresponding to the U's)
    pi <- rbeta(1,(a+sum(rn.gen)),(b+V-sum(rn.gen)))

    
    #### 9. Update xi
    which.rn.gen = which(rn.gen==1)

    
    
    len.nonzero.ui <- length(which.rn.gen)
    
    
    U.U  <- matrix(0,L*R,L*R)

    
    for(j in 1:V){ 
      
      if (j %in% which.rn.gen) {
        
        if (normalize_M) {
          
          
          U.U  <- U.U + tcrossprod(ulist[j,] / sqrt(sum(ulist[j, ] ^ 2))) 
          
        }
        
        else {
          
          U.U  <- U.U + tcrossprod(ulist[j,]) 
          
        }
        
        
      }
  
    }
    

    
    
    M  <- rinvwishart(a.wish+len.nonzero.ui, S.scale+U.U)
    
    

    
    ## update kappa1, kappa2

    U.list <- cbind(ulist1,ulist2)
    kappa.list <- rbind(kappa1,kappa2)
    Xmat.list  <- cbind(Xmat1,Xmat2)
    pi.kap.list <- rbind(pi1.kap, pi2.kap)
    
    for (kk in 1:2){
      
      Xmat.list.kk  <- Xmat.list[,((Q*(kk-1)+1):(Q*kk))]
      Xmat.list.not.kk  <- Xmat.list[,-((Q*(kk-1)+1):(Q*kk))]
      Xmat.list.not.kk1 <- Xmat.list.not.kk[,1:Q]
      
      U.list.kk     <- U.list[,(R*(kk-1)+1):(R*kk)]
      U.list.not.kk <- U.list[,-((R*(kk-1)+1):(R*kk))]
      
      kappa.list.kk <- kappa.list[kk,]
      kappa.list.not.kk <- kappa.list[-kk,]
      
      coef.temp1 <- upperTriangle(U.list.not.kk[,1:R]%*%(t(U.list.not.kk[,1:R]) * kappa.list.not.kk),
                                  diag=FALSE,byrow=T)
      
      
      
      pi.kap.kk <- pi.kap.list[(R*(kk-1)+1):(R*kk), ]
      
      for(r in 1:R){
        kappa.temp1 <- kappa.list.kk
        kappa.temp2 <- kappa.list.kk
        kappa.temp3 <- kappa.list.kk
        
        kappa.temp1[r] <- 1
        kappa.temp2[r] <- 0
        kappa.temp3[r] <- -1
        
        Byth1 <- diag(kappa.temp1)
        Byth2 <- diag(kappa.temp2)
        Byth3 <- diag(kappa.temp3)
        
        ## Upper Triangular portion of U-prime-U
        UpU.up1 <- upperTriangle(U.list.kk%*%(t(U.list.kk) * kappa.temp1),diag=FALSE,byrow=T)
        UpU.up2 <- upperTriangle(U.list.kk%*%(t(U.list.kk) * kappa.temp2),diag=FALSE,byrow=T)
        UpU.up3 <- upperTriangle(U.list.kk%*%(t(U.list.kk) * kappa.temp3),diag=FALSE,byrow=T)
        

        k.tilde <- Xmat.list.not.kk1%*%coef.temp1 + W%*%gamma +  mu
        mean.f1 <- Xmat.list.kk%*%UpU.up1
        mean.f2 <- Xmat.list.kk%*%UpU.up2
        mean.f3 <- Xmat.list.kk%*%UpU.up3
        
        prob.up1 <- sum(dnorm(k-k.tilde,mean.f1,big_omega,log=TRUE))
        prob.up2 <- sum(dnorm(k-k.tilde,mean.f2,big_omega,log=TRUE))
        prob.up3 <- sum(dnorm(k-k.tilde,mean.f3,big_omega,log=TRUE))
        prob.f.up <- c(prob.up1,prob.up2,prob.up3)
        
        #print(prob.f.up)
        if(all(prob.f.up == -Inf)) {
          kappa.list.kk[r] <- sample(c(1,0,-1),1,prob=pi.kap.kk[r,],replace=T)        
        }
        
        else {
          ind.prob <- which(prob.f.up==max(prob.f.up))[1]
          puplo <- exp(prob.f.up-prob.f.up[ind.prob])
          
          
          
          # Need to use probabilities for each different layer 

          
          kappa.list.kk[r] <- sample(c(1,0,-1),1,prob=pi.kap.kk[r,]*puplo,replace=T)
        }
        
      }
      
      kappa.list[kk, ] = kappa.list.kk
      
    }
    
    kappa1 <- kappa.list[1,]
    kappa2 <- kappa.list[2,]
    
    Byatha1 <- diag(kappa1)
    Byatha2 <- diag(kappa2)
    
    
    #### 10. Update pi.kap
    pi.kap.list <- array(NA,dim=c(L,R,3))

    for(jj in 1:L){

      
      for(r in 1:R){
        if(kappa.list[jj,r]==1){
          mm.new <- c(1,0,0)
        }else if(kappa.list[jj,r]==0){
          mm.new <- c(0,1,0)
        }else{mm.new <- c(0,0,1)}
        
        pi.kap.list[jj,r,] <- rdirichlet(1,c(a.kap[r],b.kap[r],c.kap[r])+mm.new)
      }
    }

    pi1.kap <- pi.kap.list[1,,]
    pi2.kap <- pi.kap.list[2,,]
    
    
    # Update omega
    for (nn in 1:N) {
      
      omega[nn] = rpg(num = 1, h = 1, z = mu + as.numeric(crossprod(W[nn, ], gamma)) +
                        as.numeric(crossprod(Xmat1[nn, ], coef1)) +
                        as.numeric(crossprod(Xmat2[nn, ], coef2)))
      
    }
    

    
    #### Store parameters

    mu.store[i]  <- mu
    gamma.store[i, ] <- gamma
    pi.store[i] <- pi

    
    xi.store[[i]]  <- M

    
    rn.gen.store[i,] <- rn.gen

    
    ulist.store[[i]]  <- ulist
    ulist1.store[[i]] <- ulist1
    ulist2.store[[i]] <- ulist2
    

    kappa1.store[i,]  <- kappa1
    kappa2.store[i,]  <- kappa2
    

    pi1.kap.store[i,,]  <- pi1.kap
    pi2.kap.store[i,,]  <- pi2.kap
    
    omega.store[i, ] <- omega
  }  
  
  n_samples = niter
  
  samples = seq(n_burnin+1, n_samples, by = nthinning)
  
  # Looking at posterior samples (post convergence) for each cell of B^(1), B^(2), B^(3) and B^(4)
  
  Bcells_1 = matrix(0, nrow = length(samples), ncol = ncol(Xmat1))
  Bcells_2 = matrix(0, nrow = length(samples), ncol = ncol(Xmat2))
  

  for (i in 1:length(samples)) {
    
    
    Bcells_1[i, ] = upperTriangle(ulist1.store[[samples[i]]] %*% diag(kappa1.store[samples[i], ])
                                  %*% t(ulist1.store[[samples[i]]]),
                                  diag = FALSE, byrow = T)
    
    Bcells_2[i, ] = upperTriangle(ulist2.store[[samples[i]]] %*% diag(kappa2.store[samples[i], ])
                                  %*% t(ulist2.store[[samples[i]]]),
                                  diag = FALSE, byrow = T)
    
  }
  
  coef1_MSE = NULL
  coef2_MSE = NULL
  
  if (!is.null(b1.tr) & !is.null(b2.tr)) {
  
    mean_Bcells_1 = colMeans(Bcells_1)
    mean_Bcells_2 = colMeans(Bcells_2)
    
    
    # Getting MSE for coefficients
    
    B1_coef_tr = b1.tr
    B2_coef_tr = b2.tr
    
    B1_MSE = sum((B1_coef_tr - mean_Bcells_1)^2) / sum(B1_coef_tr^2)
    
    
    B2_MSE = sum((B2_coef_tr - mean_Bcells_2)^2) / sum(B2_coef_tr^2)
    
    
    coef1_MSE = B1_MSE
    coef2_MSE = B2_MSE
    
    
    
    # Finding average coverage and length of 95% credible intervals
    
    interval_lengths1 = rep(NA, length(mean_Bcells_1))
    interval_coverages1 = rep(NA, length(mean_Bcells_1))
    
    interval_lengths2 = rep(NA, length(mean_Bcells_2))
    interval_coverages2 = rep(NA, length(mean_Bcells_2))
    
    
    for (i in 1:length(B1_coef_tr)){
      
      interval1 = quantile(Bcells_1[, i], probs = c(0.025, 0.975))
      interval_lengths1[i] = interval1[[2]] - interval1[[1]]
      interval_coverages1[i] = (interval1[[1]] < B1_coef_tr[i] &
                                  B1_coef_tr[i] < interval1[[2]])
      
      interval2 = quantile(Bcells_2[, i], probs = c(0.025, 0.975))
      interval_lengths2[i] = interval2[[2]] - interval2[[1]]
      interval_coverages2[i] = (interval2[[1]] < B2_coef_tr[i] &
                                  B2_coef_tr[i] < interval2[[2]])
    }
    
    mean_coverage1 = mean(interval_coverages1)
    mean_coverage2 = mean(interval_coverages2)
    
    mean_length1 = mean(interval_lengths1)
    mean_length2 = mean(interval_lengths2)
  }
    
  # Posterior of predictor coefficients
  
  gamma_posterior = matrix(NA, nrow = n_samples, ncol = p)
  
  for (i in 1:n_samples) {
    
    gamma_posterior[i, ] = gamma.store[samples[i], ]
    
  }
  
  
  
  # Out of Sample Prediction
  N_test = length(y_test)
  
  y_test_preds = matrix(NA, nrow = length(samples), ncol = N_test)
  
  for (i in 1:length(samples)) {
    
    y_test_preds[i, ] = exp(mu.store[samples[i]] + W_test %*% gamma_posterior[i, ] +
                            + Xmat1_test%*%Bcells_1[i, ] + Xmat2_test%*%Bcells_2[i, ]) / (1 + 
                            exp(mu.store[samples[i]] + W_test %*% gamma_posterior[i, ]
                                + Xmat1_test%*%Bcells_1[i, ] + Xmat2_test%*%Bcells_2[i, ]))
    
  }
  
  
  y_test_preds_means = colMeans(y_test_preds)
  
  auc = roc(as.numeric(y_test), y_test_preds_means)$auc
  
  R_eff1 = mean(rowSums(abs(kappa1.store[samples, ])))
  R_eff2 = mean(rowSums(abs(kappa2.store[samples, ])))
  
  node_probabilities = colMeans(rn.gen.store[samples, ])
  
  outputs = list(coef1_MSE = coef1_MSE, coef2_MSE = coef2_MSE, auc = auc, R_eff1 = R_eff1,
                 R_eff2 = R_eff2, node_probabilities = node_probabilities)
  
  return(outputs)
  
}  
  




