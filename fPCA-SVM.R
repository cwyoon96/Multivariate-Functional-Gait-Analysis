library(fda)

fpc.score = function(mat, endtime, n_pc) {
  
  pcs_df = matrix(nrow = ncol(mat), ncol = n_pc)
  
  timebasis = create.fourier.basis(c(0,endtime),as.integer(endtime * 0.8))
  
  walkfd = smooth.basis(1:endtime,mat,timebasis)
  
  walkpca = pca.fd(walkfd$fd,nharm=n_pc)
  
  col.name = c()
  
  for (i in c(1:n_pc)) {
    pcs_df[,i] = walkpca$scores[,i]
    col.name = c(col.name, paste('pcs_',i, sep = ''))
  }
  print(walkpca$varprop)
  print(sum(walkpca$varprop))

  par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
  plot(walkpca$values[1:8],xlab='component',ylab='variance',col="red",
       cex.lab=1.5,cex.axis=1.5,cex=2)
  
  colnames(pcs_df) = col.name
  
  return(data.frame(pcs_df))
}

SVM_CV = function(data, n_cv, manual = FALSE, manual_index = NaN, positive = '1'){
  
  if(manual == TRUE){
    train = data[manual_index,]
    test = data[-manual_index,]
  } else {
    trainIndex = createDataPartition(as.factor(data$gmfcs), p = .7, list = FALSE)
    train = data[trainIndex,]
    test = data[-trainIndex,]
  }
  
  cvIndex = createFolds(as.factor(train$gmfcs), k = n_cv, list = FALSE)
  train['cv'] = cvIndex
  
  param_list = matrix(nrow = 20, ncol = 2)
  param_list[,1] = seq(0.1,5,length.out = 20)
  
  for (j in seq(1,20)){
    err_mean = 0
    for (i in c(1:n_cv)) {
      train_cv = train[train$cv != i, ]
      test_cv = train[train$cv ==i, ]
      
      md = svm(formula = gmfcs ~ pcs_1 + pcs_2 + pcs_3, data = train_cv, type = "C-classification",kernel = 'linear', cost = param_list[j,1])
      cm = confusionMatrix(reference = test_cv[,'gmfcs'],
                           data = predict(md,test_cv[,c('pcs_1','pcs_2','pcs_3')]), 
                           positive = positive,mode = "everything")
      tab = cm$table
      err = sum(tab) - sum(diag(tab))
      
      err_mean = err_mean + err/n_cv
    }
    
    param_list[j,2] = err_mean
    
  }
  
  opt_param = which(param_list[,2] == min(param_list[,2]))
  
  md = svm(formula = gmfcs ~ pcs_1 + pcs_2 + pcs_3, data = train, type = "C-classification",kernel = 'linear', cost = param_list[opt_param[1],1])
  
  cm = confusionMatrix(reference = test[,'gmfcs'],
                       data = predict(md,test[,c('pcs_1','pcs_2','pcs_3')]), 
                       positive = positive, mode = "everything")
  
  return(cm)
  
}
