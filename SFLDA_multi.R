library(pracma)
library(VGAM)
library(MASS)

is_constant = function(xb){
  const = c(var(xb[,1]), var(xb[,2]))
  if (sum(const) == 0){
    return(NULL)
  } else{
    return(xb[,which(const != 0),drop = FALSE])
  }
}

SDelta_multi = function(xtr1,xtr2,xtr3) {
  # Compute pooled S
  ntr1 = nrow(xtr1)
  ntr2 = nrow(xtr2)
  ntr3 = nrow(xtr3)
  p = ncol(xtr1)
  
  S = (ntr1 - 1)*cov(xtr1) + (ntr2 - 1)*cov(xtr2) + (ntr3 - 1)*cov(xtr3) 
  S = S/(ntr1 + ntr2 + ntr3 - 3)
  
  # Compute X_bar
  mu1 = as.matrix(colMeans(xtr1))
  mu2 = as.matrix(colMeans(xtr2))
  mu3 = as.matrix(colMeans(xtr3))
  
  # Compute cov(t(A)) (pxp matrix)
  
  A = cbind(mu1,mu2,mu3)
  A_cov = cov(t(A))
  
  # Eigenvalue Decomposition and pick first 2 eigenvectors as deltas
  
  eigen_A_cov = eigen(A_cov)
  
  delta1 = eigen_A_cov$vectors[,1]
  delta2 = eigen_A_cov$vectors[,2]
  
  return(list('S' = S, 'delta1' = delta1, 'delta2' = delta2))
}

SFLDA_multi_beta = function(S, delta1, delta2,tau, lambda, orthogonal = TRUE, multivariate = FALSE){
  
  # Compute each beta with its delta
  
  p = length(delta1)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))
  
  if (multivariate == TRUE){
    D[100,] = rep(0,300)
    D[200,] = rep(0,300)
  }
  
  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')
  
  S = S / norm(S, type = 'F')
  
  b1 = coordlasso((1-tau)*S + tau*DD, delta1, lambda, 1e-8, 500)
  b2 = coordlasso((1-tau)*S + tau*DD, delta2, lambda, 1e-8, 500)
  
  # make each beta orthogonal
  
  beta = cbind(b1,b2)
  
  if ((pracma::Rank(beta) == 2) & (orthogonal == TRUE)){
    beta = gramSchmidt(beta)$Q
  }
  
  
  return(beta)
}

SFLDA_multi = function(xtr1, xtr2, xtr3, xte1,xte2,xte3, tau, lambda, multivariate = FALSE){
  
  ntr1 = nrow(xtr1)
  ntr2 = nrow(xtr2)
  ntr3 = nrow(xtr3)
  
  nte1 = nrow(xte1)
  nte2 = nrow(xte2)
  nte3 = nrow(xte3)
  
  SD = SDelta_multi(xtr1,xtr2,xtr3)
  
  S = SD$S
  
  delta1 = SD$delta1
  delta2 = SD$delta2
  
  beta = SFLDA_multi_beta(S, delta1, delta2, tau, lambda, orthogonal = TRUE, multivariate = multivariate)
  
  xb = rbind(xtr1%*%beta, xtr2%*%beta, xtr3%*%beta)

  xb = data.frame(xb)
  colnames(xb) = c('xb1','xb2')
  
  if (is.null(is_constant(xb))){
    class = rep(2,(nte1 + nte2 + nte3))
    class = as.factor(class)
    
    xb_te = rbind(xte1%*%beta, xte2%*%beta, xte3%*%beta)
    xb_te = data.frame(xb_te)
    
    colnames(xb_te) = c('xb1','xb2')
    
    gmfcs_te = c(rep(1,nte1),rep(2,nte2),rep(3,nte3))
    xb_te$gmfcs = as.factor(gmfcs_te)
    
  } else {
    xb = is_constant(xb)
    gmfcs = c(rep(1,ntr1),rep(2,ntr2),rep(3,ntr3))
    xb$gmfcs = as.factor(gmfcs)
    
    md = lda(formula = gmfcs ~ ., data = xb) 
    
    xb_te = rbind(xte1%*%beta, xte2%*%beta, xte3%*%beta)
    xb_te = data.frame(xb_te)
    colnames(xb_te) = c('xb1','xb2')
    
    xb_te = is_constant(xb_te)
    gmfcs_te = c(rep(1,nte1),rep(2,nte2),rep(3,nte3))
    xb_te$gmfcs = as.factor(gmfcs_te)
    
    pred = predict(md, within(xb_te, rm(gmfcs)))
    
    class = pred$class
    
  }
  
  return(list('class' = class, 'yte' = xb_te$gmfcs, 'beta' = beta, 'xb_tr' = xb, 'xb_te' = xb_te))
  
}
