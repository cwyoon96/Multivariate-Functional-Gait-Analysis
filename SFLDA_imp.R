clean_mat = function(df) {
  mat = as.matrix(df[1:(length(df)-2)])
  colnames(mat) = NULL
  rownames(mat) = NULL
  return(mat)
}

clean_mat_test = function(df) {
  mat = as.matrix(df[1:(length(df)-1)])
  colnames(mat) = NULL
  rownames(mat) = NULL
  return(mat)
}

SFLDA_CV = function(data, n_cv, manual = FALSE, manual_index = NaN) {
  if(manual == TRUE){
    train = data[manual_index,]
    test = data[-manual_index,]
  } else {
    trainIndex = createDataPartition(as.factor(data$group), p = .7, list = FALSE)
    train = data[trainIndex,]
    test = data[-trainIndex,]
  }
  
  cvIndex = createFolds(as.factor(train$group), k = n_cv, list = FALSE)
  train['cv'] = cvIndex
  
  params = matrix(nrow = 100,ncol = 3)
  
  tau = 2^seq(-4,-0.1,length.out = 10)
  lambda = 2^seq(-6,-1, length.out = 10)
  
  params[,1] = rep(tau,10)
  params[,2] = rep(lambda, each = 10)
  
  pb = progress_bar$new(total = 100)
  
  for (j in c(1:100)){
    pb$tick()
    
    err_mean = 0
    
    for (i in c(1:n_cv)) {
      train_cv = train[train$cv != i, ]
      test_cv = train[train$cv ==i, ]
      
      train_cv0 = clean_mat(train_cv[train_cv$group == 0, ])
      train_cv1 = clean_mat(train_cv[train_cv$group == 1, ])
      test_cv0 = clean_mat(test_cv[test_cv$group == 0, ])
      test_cv1 = clean_mat(test_cv[test_cv$group ==1, ])
      
      trcv.SDelta = SDelta(train_cv0, train_cv1)
      tecv.SDelta = SDelta(test_cv0, test_cv1)
      
      
      SFLDAcv.value = SFLDA_binary(test_cv0,test_cv1,trcv.SDelta$S,trcv.SDelta$mu0,trcv.SDelta$mu1,params[j,1],params[j,2],
                                   tecv.SDelta$S,tecv.SDelta$mu0,tecv.SDelta$mu1)
      
      err_mean = err_mean + SFLDAcv.value$err/n_cv
      
    }
    params[j,3] = err_mean
  }
  
  params = data.frame(params)
  colnames(params) = c('tau.','lambda.','err.')
  attach(params)
  params = params[order(err.,-lambda.,-tau.),]
  detach(params)
  
  train0 = clean_mat(train[train$group == 0,])
  train1 = clean_mat(train[train$group == 1,])
  test0 = clean_mat_test(test[test$group == 0,])
  test1 = clean_mat_test(test[test$group == 1,])
  
  tr.SDelta = SDelta(train0, train1)
  te.SDelta = SDelta(test0, test1)
  
  SFLDA.value = SFLDA_binary(test0,test1,tr.SDelta$S,tr.SDelta$mu0,tr.SDelta$mu1,params[1,1],params[1,2],
                             te.SDelta$S,te.SDelta$mu0,te.SDelta$mu1)
  
  result = list('b' = SFLDA.value$b, 'acc' = SFLDA.value$acc,
                'class' = SFLDA.value$class, 'yte' = SFLDA.value$yte,
                'tr0' = train0%*%(SFLDA.value$b), 'tr1' = train1%*%(SFLDA.value$b),
                'te0' = test0%*%(SFLDA.value$b), 'te1' = test1%*%(SFLDA.value$b),
                'pm0' = t(tr.SDelta$mu0)%*%(SFLDA.value$b), 'pm1' = t(tr.SDelta$mu1)%*%(SFLDA.value$b),
                'tau' = params[1,1], 'lambda' = params[1,2])

  return(result)
  
}

SFLDA_multi_CV = function(data, n_cv, manual = FALSE, manual_index = NaN, multivariate = FALSE) {
  if(manual == TRUE){
    train = data[manual_index,]
    test = data[-manual_index,]
  } else {
    trainIndex = createDataPartition(as.factor(data$group), p = .8, list = FALSE)
    train = data[trainIndex,]
    test = data[-trainIndex,]
  }
  
  cvIndex = createFolds(as.factor(train$group), k = n_cv, list = FALSE)
  train['cv'] = cvIndex
  
  params = matrix(nrow = 100,ncol = 3)
  
  tau = 2^seq(-4,-0.1,length.out = 10)
  lambda = 2^seq(-6,-1, length.out = 10)
  
  params[,1] = rep(tau,10)
  params[,2] = rep(lambda, each = 10)
  
  pb = progress_bar$new(total = 100)
  
  for (j in c(1:100)){
    pb$tick()
    
    err_mean = 0
    
    for (i in c(1:n_cv)) {
      train_cv = train[train$cv != i, ]
      test_cv = train[train$cv ==i, ]
      
      train_cv1 = clean_mat(train_cv[train_cv$group == 1, ])
      train_cv2 = clean_mat(train_cv[train_cv$group == 2, ])
      train_cv3 = clean_mat(train_cv[train_cv$group == 3, ])
      
      test_cv1 = clean_mat(test_cv[test_cv$group == 1, ])
      test_cv2 = clean_mat(test_cv[test_cv$group == 2, ])
      test_cv3 = clean_mat(test_cv[test_cv$group == 3, ])
      
      
      SFLDAcv_multi.value = SFLDA_multi(train_cv1, train_cv2, train_cv3, test_cv1, test_cv2, test_cv3,
                                        tau = params[j,1], lambda =  params[j,2], multivariate = multivariate)
      class = SFLDAcv_multi.value$class
      yte = SFLDAcv_multi.value$yte
      
      cm = confusionMatrix(data = class, reference = yte, positive = '3')
      
      tab = cm$table
      
      err = sum(tab) - sum(diag(tab))
      
      err_mean = err_mean + err/n_cv
      
    }
    params[j,3] = err_mean
  }
  
  params = data.frame(params)
  colnames(params) = c('tau.','lambda.','err.')
  params = params[order(params[,'err.'],-params[,'lambda.'],-params[,'tau.']),]
  
  train1 = clean_mat(train[train$group == 1,])
  train2 = clean_mat(train[train$group == 2,])
  train3 = clean_mat(train[train$group == 3,])
  
  test1 = clean_mat_test(test[test$group == 1,])
  test2 = clean_mat_test(test[test$group == 2,])
  test3 = clean_mat_test(test[test$group == 3,])
  
  SFLDA_multi.value = SFLDA_multi(train1, train2, train3, test1, test2, test3,
                                  tau = params[1,1],
                                  lambda = params[1,2],
                                  multivariate = multivariate)
  
  result = list('beta' = SFLDA_multi.value$beta, 'class' = SFLDA_multi.value$class, 'yte' = SFLDA_multi.value$yte,
                'xb_tr' = SFLDA_multi.value$xb_tr, 'xb_te' = SFLDA_multi.value$xb_te)
  
  return(result)
  
}

