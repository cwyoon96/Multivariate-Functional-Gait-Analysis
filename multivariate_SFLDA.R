### Methods Section
setwd("/Volumes/GoogleDrive/My Drive/KAIST/HS Lab/Projects/Ajou Hospital")

library(progress)
library(ggplot2)
library(dplyr)
library(e1071)
library(caret)
library(reshape2)

source('./Paper/R codes/common_func.R')
source('./Paper/R codes/fPCA-SVM.R')
source('./Paper/R codes/flr.R')
source('./Paper/R codes/SFLDA.R')
source('./Paper/R codes/SFLDA_multi.R')
source('./Paper/R codes/SFLDA_imp.R')

data = read.csv("./Paper/data/preprocessed_data_2.csv")
columns = colnames(data)[5:16]

### normal vs patient
data_patient = data

data_patient[data_patient$gmfcs >= 1,'gmfcs'] = 1
rownames(data_patient) = NULL

### among each level

pick_two = function(data, a,b){ # pick a and b gmfcs level
  df = data[data$gmfcs %in% c(a,b),]
  df[df$gmfcs == a,'gmfcs'] = 0
  df[df$gmfcs == b,'gmfcs'] = 1
  return(df)
}

# 0 vs 1

data_01 = pick_two(data,0,1)
rownames(data_01) = NULL

# 1 vs 2

data_12 = pick_two(data,1,2)
rownames(data_12) = NULL

# 2 vs 3

data_23 = pick_two(data,2,3)
rownames(data_23) = NULL

multivariate_SFLDA_binary = function(xte0, xte1, S, mu0, mu1, tau, lambda){
  
  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))
  
  D[100,] = rep(0,300)
  D[200,] = rep(0,300)
  
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
  
  return(list('b' = b, 'err' = err, 'acc' = acc, 'class' = class, 'yte' = yte))
}

multivariate_SFLDA_get_b = function(S, mu0, mu1, tau, lambda){
  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))
  
  D[100,] = rep(0,300)
  D[200,] = rep(0,300)
  
  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')
  
  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')
  
  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)
  
  return(b)
}

multivariate_SFLDA_CV = function(data, n_cv, manual = FALSE, manual_index = NaN) {
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
      
      
      SFLDAcv.value = multivariate_SFLDA_binary(test_cv0,test_cv1,trcv.SDelta$S,trcv.SDelta$mu0,trcv.SDelta$mu1,params[j,1],params[j,2])
                                   
      
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
  
  SFLDA.value = multivariate_SFLDA_binary(test0,test1,tr.SDelta$S,tr.SDelta$mu0,tr.SDelta$mu1,params[1,1],params[1,2])
  
  result = list('b' = SFLDA.value$b, 'acc' = SFLDA.value$acc,
                'class' = SFLDA.value$class, 'yte' = SFLDA.value$yte,
                'tr0' = train0%*%(SFLDA.value$b), 'tr1' = train1%*%(SFLDA.value$b),
                'te0' = test0%*%(SFLDA.value$b), 'te1' = test1%*%(SFLDA.value$b),
                'pm0' = t(tr.SDelta$mu0)%*%(SFLDA.value$b), 'pm1' = t(tr.SDelta$mu1)%*%(SFLDA.value$b),
                'tau' = params[1,1], 'lambda' = params[1,2])
  
  return(result)
  
}

perf_matrix = function(){
  perf = vector('list',2)
  names(perf) = c('R','L')
  
  for (name in names(perf)){
    perf[[name]] = vector('list',10)
  }
  
  return(perf)
}



### SFLDA

multivariate_SFLDA_comparision = function(data){
  SFLDA_perf = perf_matrix()
  
  gmfcs = data$gmfcs[seq(1,nrow(data)-99,100)]
  
  set.seed(123)
  partition_index = createDataPartition(as.factor(gmfcs), p = .7, 
                                        list = FALSE, 
                                        times = 10)
  
  beta_list = vector('list',2)
  names(beta_list) = c('R', 'L')
  
  param_list = vector('list',2)
  names(param_list) = c('R', 'L')
  
  for (direc in c('R','L')){
    beta_list[[direc]] = matrix(nrow = 300, ncol = 10)
    param_list[[direc]] = matrix(nrow = 10, ncol = 2)
    colnames(param_list[[direc]]) = c('tau','lambda')
  }
  
  for (direc in c('R','L')){
    
    hip_data = t(to_matrix(data, paste(direc,'HIP_Flex_ANG',sep ='_')))
    knee_data = t(to_matrix(data,paste(direc,'KNEE_Flex_ANG',sep ='_')))
    ankle_data = t(to_matrix(data,paste(direc,'ANK_Flex_ANG',sep ='_')))
    
    df = cbind(hip_data, knee_data, ankle_data)
    
    df = data.frame(df)
    
    df['group'] = gmfcs
    
    for (i in seq(1:10)){
      result = multivariate_SFLDA_CV(df, 5, manual = TRUE, manual_index = partition_index[,i])
      print(paste(i,'th cv done.', sep = ''))
      
      beta_list[[direc]][,i] = result$b
      
      param_list[[direc]][i,] = c(result$tau, result$lambda)
      
      cm = confusionMatrix(reference = as.factor(result$yte), data = as.factor(result$class), positive = '2')
      
      SFLDA_perf[[direc]][[i]] = cm
      
    }
    print(paste(direc, 'done'))
  }
  
  return(list('SFLDA' = SFLDA_perf, 'b_list' = beta_list, 'p_list' = param_list))
}

multivariate_SFLDA_01 = multivariate_SFLDA_comparision(data_01)
multivariate_SFLDA_12 = multivariate_SFLDA_comparision(data_12)
multivariate_SFLDA_23 = multivariate_SFLDA_comparision(data_23)
multivariate_SFLDA_patient = multivariate_SFLDA_comparision(data_patient)



calculate_metric = function(lst){
    
  metric = matrix(nrow = 6, ncol = 2)
  colnames(metric) = c('R', 'L')
  rownames(metric) = c('acc','acc_std','fnr','fnr_std','for','for_std')
    
  for (direc in c('R','L')){
    acc = c()
    fnr = c()
    for_ = c()
    for (i in seq(1,10)){
      acc = c(acc, lst$SFLDA[[direc]][[i]]$overall['Accuracy'])
      fnr = c(fnr, 1 - lst$SFLDA[[direc]][[i]]$byClass['Sensitivity'])
      for_ = c(for_, 1 - lst$SFLDA[[direc]][[i]]$byClass['Neg Pred Value'])
    }
      
    metric['acc',direc] = mean(acc,na.rm = TRUE)
    metric['acc_std',direc] = sd(acc, na.rm = TRUE)
      
    metric['fnr',direc] = mean(fnr, na.rm = TRUE)
    metric['fnr_std',direc] = sd(fnr, na.rm = TRUE)
      
    metric['for',direc] = mean(for_, na.rm = TRUE)
    metric['for_std',direc] = sd(for_, na.rm = TRUE)
    
  }

  
  return(metric)
}

round(calculate_metric(multivariate_SFLDA_01),2)
round(calculate_metric(multivariate_SFLDA_12),2)
round(calculate_metric(multivariate_SFLDA_23),2)
round(calculate_metric(multivariate_SFLDA_patient),2)


save(multivariate_SFLDA_01, file = "./Paper/data/multivariate_SFLDA_01.RData")
save(multivariate_SFLDA_12, file = "./Paper/data/multivariate_SFLDA_12.RData")
save(multivariate_SFLDA_23, file = "./Paper/data/multivariate_SFLDA_23.RData")
save(multivariate_SFLDA_patient, file = "./Paper/data/multivariate_SFLDA_patient.RData")

get_best_param = function(lst, direc){
  acc = c()
  for (i in seq(1:10)){
    acc = c(acc,lst$SFLDA[[direc]][[i]]$overall[1])
  }
  max_cv = which(acc == max(acc))
  
  return(matrix(lst$p_list[[direc]][max_cv,],ncol = 2))
}


get_b_mat = function(data,lst){
  b_matrix = matrix(nrow = 300, ncol = 2)
  colnames(b_matrix) = c('R','L')
  for (direc in c('R','L')){
    
    total_0 = cbind(t(to_matrix(data[data$gmfcs == 0,],paste(direc,'HIP_Flex_ANG',sep = '_'))), 
                  t(to_matrix(data[data$gmfcs == 0,],paste(direc,'KNEE_Flex_ANG',sep = '_'))),
                  t(to_matrix(data[data$gmfcs == 0,],paste(direc,'ANK_Flex_ANG',sep = '_'))))
    
    total_1 = cbind(t(to_matrix(data[data$gmfcs == 1,],paste(direc,'HIP_Flex_ANG',sep = '_'))), 
                    t(to_matrix(data[data$gmfcs == 1,],paste(direc,'KNEE_Flex_ANG',sep = '_'))),
                    t(to_matrix(data[data$gmfcs == 1,],paste(direc,'ANK_Flex_ANG',sep = '_'))))
    
    sdelta = SDelta(total_0, total_1)
    param = get_best_param(lst, direc)[1,]
    b_matrix[,direc] = multivariate_SFLDA_get_b(sdelta$S,sdelta$mu0,sdelta$mu1, param[1], param[2])
  }
  
  return(b_matrix)
}

multivariate_b_matrix_01 = get_b_mat(data_01, multivariate_SFLDA_01)
multivariate_b_matrix_12 = get_b_mat(data_12, multivariate_SFLDA_12)
multivariate_b_matrix_23 = get_b_mat(data_23, multivariate_SFLDA_23)
multivariate_b_matrix_patient = get_b_mat(data_patient, multivariate_SFLDA_patient)

combine_multivarite_plot = function(mat1,mat2,mat3, direc){
  hip = cbind(mat1[1:100,direc],mat2[1:100,direc],mat3[1:100,direc])
  knee = cbind(mat1[101:200,direc],mat2[101:200,direc],mat3[101:200,direc])
  ankle = cbind(mat1[201:300,direc],mat2[201:300,direc],mat3[201:300,direc])
  
  return(list(hip = hip, knee = knee, ankle = ankle))
}

Right_part = combine_multivarite_plot(multivariate_b_matrix_01,multivariate_b_matrix_12,multivariate_b_matrix_23, 'R')
Left_part = combine_multivarite_plot(multivariate_b_matrix_01,multivariate_b_matrix_12,multivariate_b_matrix_23, 'L')

for (loc in c('hip','knee','ankle')){
  b_mat  = data.frame(cbind(Right_part[[loc]],seq(1,100)))
  colnames(b_mat) = c('0 vs. 1','1 vs. 2','2 vs. 3', 'time')
  
  p = ggplot(data = melt(b_mat, id = 'time'), aes(x = time, y = value, group = variable, color = variable)) +
    geom_line(size = 3, aes(linetype = variable)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'None',
          text = element_text(size = 25),
          panel.background = element_rect(colour = "black", size=2)) +
    scale_x_continuous(limits=c(0, 101), expand = c(0, 0)) +
    labs(y = expression(hat(beta)),
         x = 'Gait cycle') +
    scale_color_brewer(palette="Accent") +
    scale_linetype_manual(values=c("solid","longdash", "twodash")) +
    scale_y_continuous(limits = c(-0.6, 0.4))
  
  print(p)
  ggsave(filename =paste('./Paper/Plots/Methods/multivariate_SFLDA_beta_R_',loc,'.png',sep = ''), plot = p,width = 10.47, height = 6.98)
  
}

for (loc in c('hip','knee','ankle')){
  b_mat  = data.frame(cbind(Left_part[[loc]],seq(1,100)))
  colnames(b_mat) = c('0 vs. 1','1 vs. 2','2 vs. 3', 'time')
  
  p = ggplot(data = melt(b_mat, id = 'time'), aes(x = time, y = value, group = variable, color = variable)) +
    geom_line(size = 3, aes(linetype = variable)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'None',
          text = element_text(size = 25),
          panel.background = element_rect(colour = "black", size=2)) +
    scale_x_continuous(limits=c(0, 101), expand = c(0, 0)) +
    labs(y = expression(hat(beta)),
         x = 'Gait cycle') +
    scale_color_brewer(palette="Accent") +
    scale_linetype_manual(values=c("solid","longdash", "twodash")) +
    scale_y_continuous(limits = c(-0.6, 0.4))
  
  print(p)
  ggsave(filename =paste('./Paper/Plots/Methods/multivariate_SFLDA_beta_L_',loc,'.png',sep = ''), plot = p,width = 10.47, height = 6.98)
  
  
}

R_hip = multivariate_b_matrix_patient[1:100,'R']
R_knee = multivariate_b_matrix_patient[101:200,'R']
R_ankle = multivariate_b_matrix_patient[201:300,'R']

L_hip = multivariate_b_matrix_patient[1:100,'L']
L_knee = multivariate_b_matrix_patient[101:200,'L']
L_ankle = multivariate_b_matrix_patient[201:300,'L']

multivariate_patient = list(R_hip, R_knee, R_ankle, L_hip, L_knee, L_ankle)

for (i in seq(1,6)){
  b_mat  = data.frame(cbind(multivariate_patient[[i]], seq(1,100)))
  colnames(b_mat) = c('0 vs. {1,2,3}', 'time')
  
  g = ggplot(data = melt(b_mat, id = 'time'), aes(x = time, y = value, group = variable, color = variable)) +
    geom_line(size = 3) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'None',
          text = element_text(size = 25),
          panel.background = element_rect(colour = "black", size=2)) +
    scale_x_continuous(limits=c(0, 101), expand = c(0, 0)) +
    labs(y = expression(hat(beta)),
         x = 'Gait cycle') +
    scale_y_continuous(limits = c(-0.6, 0.4))
  
  print(g)
  
  ggsave(filename =paste('./Paper/Plots/Methods/multivariate_SFLDA_beta_patient_',as.character(i),'.png',sep = ''), plot = g,width = 10.47, height = 6.98)
}

### Multivaraite SFLDA multiclass

# SFLDA_multi

multi_sen = function(tab){
  inv_sen.3 = (tab[1,3] + tab[2,3]) /sum(tab[,3])
  inv_sen.2 = tab[1,2]/sum(tab[,2])
  
  return(mean(c(inv_sen.3,inv_sen.2)))
}

multi_npv = function(tab){
  inv_npv.3 = (tab[1,2] + tab[1,3]) /sum(tab[1,])
  inv_npv.2 = tab[2,3]/sum(tab[2,])
  
  return(mean(c(inv_npv.3,inv_npv.2)))
}

data_multi = data[data$gmfcs > 0,]

rownames(data_multi) = NULL

gmfcs = data_multi$gmfcs[seq(1,nrow(data_multi)-99,100)]

set.seed(123)
partition_index = createDataPartition(as.factor(gmfcs), p = .7, 
                                      list = FALSE, 
                                      times = 10)

SFLDA_perf = vector('list', length = 3)
names(SFLDA_perf) = c('acc','1-sen','1-npv')

for (metric in c('acc','1-sen','1-npv')){
  SFLDA_perf[[metric]] = matrix(nrow = 10, ncol = 2)
  colnames(SFLDA_perf[[metric]]) = c('R','L')
}


for (direc in c('L')){
  
  hip_data = t(to_matrix(data_multi, paste(direc,'HIP_Flex_ANG',sep ='_')))
  knee_data = t(to_matrix(data_multi,paste(direc,'KNEE_Flex_ANG',sep ='_')))
  ankle_data = t(to_matrix(data_multi,paste(direc,'ANK_Flex_ANG',sep ='_')))
  
  df = cbind(hip_data, knee_data, ankle_data)
  
  df = data.frame(df)
  
  df['group'] = gmfcs
  
  for (i in seq(1:10)){
    result = SFLDA_multi_CV(df, 5, manual = TRUE, manual_index = partition_index[,i], multivariate = TRUE)
    print(paste(i,'th cv done.', sep = ''))
    
    cm = confusionMatrix(reference = as.factor(result$yte), data = as.factor(result$class), positive = '3')
    
    SFLDA_perf[['acc']][i,direc] = cm$overall[1]
    SFLDA_perf[['1-sen']][i,direc] = multi_sen(cm$table)
    SFLDA_perf[['1-npv']][i,direc] = multi_npv(cm$table)
    
  }
  print(paste(direc, 'done'))
}

save(SFLDA_perf, file = "./Paper/data/multivariate_SFLDA_multiclass.RData")

