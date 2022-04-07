library(fda.usc)

make_flr = function(train,y, n_basis){
  basis.1 = create.fourier.basis(c(1,100),n_basis)
  basis.2 = create.fourier.basis(c(1,100),n_basis)
  
  basis.x = list('x' = basis.1)
  basis.b = list('x' = basis.2)
  
  x = fdata(train)
  y = y
  
  ldata = list('df' = y, 'x' = x)
  
  flr=classif.glm(y~x,data=ldata,basis.x=basis.x,
                  basis.b=basis.b)
  
  return(flr)
}

FLR_CV = function(data,gmfcs, n_cv, manual = FALSE, manual_index = NaN, positive = '1'){
  
  if(manual == TRUE){
    train = data[manual_index,]
    test = data[-manual_index,]
  } else {
    trainIndex = createDataPartition(as.factor(gmfcs), p = .7, list = FALSE)
    train = data[trainIndex,]
    test = data[-trainIndex,]
  }
  
  gmfcs_train = gmfcs[manual_index]
  gmfcs_test = gmfcs[-manual_index]
  
  cvIndex = createFolds(as.factor(gmfcs_train), k = n_cv, list = FALSE)
  
  params = matrix(nrow = 10, ncol = 2)
  
  params[,1] = seq(from = 3, by = 5, length.out = 10)
  
  
  pb = progress_bar$new(total = 10)
  
  for (j in seq(1,10)){
    pb$tick()
    
    err_mean = 0
    
    for (i in c(1:n_cv)) {
      train_cv = train[cvIndex != i, ]
      test_cv = train[cvIndex == i, ]
      
      y = gmfcs_train[cvIndex != i]
      y = as.data.frame(y)
      
      flr = make_flr(train_cv,y, n_basis = params[j,1])
      
      valid = fdata(test_cv)
      
      newdata = list('x' = valid)
      
      pred = predict(flr, newdata)
      
      cm = confusionMatrix(reference = as.factor(gmfcs_train[cvIndex == i]), data = pred, positive = positive)
      
      tab = cm$table
      err = sum(tab) - sum(diag(tab))
      
      err_mean = err_mean + err/n_cv
    }
    
    params[j,2] = err_mean
    
  }
  
  opt_param = which(params[,2] == min(params[,2]))

  y = gmfcs_train
  y = as.data.frame(y)
  
  flr = make_flr(train,y, params[opt_param[1],1])
  
  coef = flr$fit[[1]]$coefficients
  
  test_set = fdata(test)
  
  newdata = list('x' = test_set)
  
  pred = predict(flr, newdata)
  
  cm = confusionMatrix(reference = as.factor(gmfcs_test), data = pred, positive = positive)
  
  return(list('cm' = cm, 'coef' = coef, 'n_basis' = params[opt_param[1],1]))
  
}

