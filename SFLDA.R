vector_norm = function(x) {
  return(sqrt(sum(x^2)))
}

pdist2 = function(x,y) {
  dist.vec = matrix(nrow = nrow(x), ncol = ncol(x))
  for (i in c(1:nrow(x))) {
    dist.vec[i] = sqrt((x[i] - y)^2)
  }
  return(dist.vec)
}

wthresh = function(x,sorh,t){
  # WTHRESH Perform soft or hard thresholding.
  
  if (sorh == 's') {
    tmp = (abs(x) - t)
    tmp = (tmp+abs(tmp))/2
    y = sign(x) * tmp
  } else if (sorh == 'h') {
    y = x * (abs(x) > t)
  } else {
    print('Invalid argument value.')
  }
  return(y)
}

coordlasso = function(A, delta, lambda, tol, maxiter) {
  # Find solution
  
  # minimize 1/2 b'Ab - delta'b + m ||b||_1
  # hat b_l = S_m(delta_l - sum_{j != l} a_{lj} b_j)/a_{ll}
  # initial value
  
  p = nrow(A)
  
  b_ini = solve(A)%*%delta
  b_ini = b_ini/vector_norm(b_ini)
  
  bnew = b_ini
  
  converged = FALSE
  step = 0
  dA = as.matrix(diag(A))
  
  # obj_new = 1/2 *bnew'*A*bnew - delta'*bnew + lambda * sum(abs(bnew));
  
  while ((step < maxiter) & !converged) {
    bold = bnew
    step = step + 1
    
    for (j in c(1:p)) {
      x = delta[j] - A[j,]%*%bold + dA[j]*bold[j]
      
      x = wthresh(x,'s',lambda)
      x = x/dA[j]
      
      bold[j] = x
    }
    conv_criterion = vector_norm(bnew - bold)/vector_norm(bold)
    converged = conv_criterion < tol
    
    bnew = bold
    
    if (vector_norm(bnew) == 0) break
    
  }
  
  b = bnew
  
  if (vector_norm(b) != 0) {
    b = b / vector_norm(b)
  }
  
  return(b)
}

centroid_class = function(V,m1,m2, Xtu, Ytu) {
  # make a class assignment given the direction vector and centers
  
  # m1, m2 : 1 x p,, means of training data
  # V: p x r, PLS directions
  # vaX, vaY: n x p, n x 1 data to classify
  
  pX = Xtu %*% V
  
  pm1 = t(m1) %*% V
  pm2 = t(m2) %*% V
  
  D1 = pdist2(pX,pm1)
  D2 = pdist2(pX,pm2)
  
  Class = (D1 > D2) + 1
  
  Err = sum(Class != Ytu)
  
  CErr = NaN
  
  return(list('Err'=Err, 'Class' = Class))
}

SDelta = function(xtr0cv,xtr1cv) {
  ntr0cv = nrow(xtr0cv)
  ntr1cv = nrow(xtr1cv)
  p = ncol(xtr1cv)
  
  S = (ntr0cv - 1)*cov(xtr0cv) + (ntr1cv - 1)*cov(xtr1cv)
  S = S/(ntr0cv + ntr1cv - 2)
  mu0 = as.matrix(colMeans(xtr0cv))
  mu1 = as.matrix(colMeans(xtr1cv))
  
  return(list('S' = S, 'mu0' = mu0, 'mu1' = mu1))
}

SFLDA_binary = function(xte0, xte1, S, mu0, mu1, tau, lambda, Ste, mu0te, mu1te){
  
  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))
  
  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')
  
  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')
  
  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)
  
  nte0 = nrow(xte0)
  nte1 = nrow(xte1)
  
  xte = rbind(xte0,xte1)
  yte = as.matrix(c(rep(1,nte0), rep(2,nte1)))
  
  cclass = centroid_class(b,mu0,mu1,xte,yte)
  err = cclass$Err
  class = cclass$Class
  
  acc = 1 - (err/nrow(xte))
  
  if (all(Ste == 0)) {
    obj = 100
  } else {
    delta_te = (mu0te - mu1te)
    
    obj = 0.5 * t(b) %*% Ste %*%b - t(b) %*% delta_te
  }
  
  return(list('b' = b, 'err' = err, 'acc' = acc, 'obj' = obj, 'class' = class, 'yte' = yte))
}

SFLDA_get_b = function(S, mu0, mu1, tau, lambda) {
  
  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))
  
  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')
  
  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')
  
  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)
  
  return(b)
}