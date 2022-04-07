### Methods Section
setwd("/Volumes/GoogleDrive/My Drive/KAIST/HS Lab/Projects/CP_SFLDA")

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


perf_matrix = function(){
  perf = vector('list',6)
  names(perf) = c('R_HIP_Flex_ANG','L_HIP_Flex_ANG','R_KNEE_Flex_ANG','L_KNEE_Flex_ANG','R_ANK_Flex_ANG','L_ANK_Flex_ANG')
  
  for (name in names(perf)){
    perf[[name]] = vector('list',10)
  }
  
  return(perf)
}

multi_perf_matrix = function(){
  perf = vector('list',3)
  names(perf) = c('acc','1-sen','1-npv')
  
  for (name in names(perf)){
    perf[[name]] = matrix(nrow = 10,ncol = 6)
    colnames(perf[[name]]) = c('R_HIP_Flex_ANG','L_HIP_Flex_ANG','R_KNEE_Flex_ANG','L_KNEE_Flex_ANG','R_ANK_Flex_ANG','L_ANK_Flex_ANG')
  }
  
  return(perf)
}


# load data
data = read.csv("./Paper/data/preprocessed_data_2.csv")
columns = colnames(data)[5:16]

model_comparsion = function(data){
  
  gmfcs = data$gmfcs[seq(1,nrow(data)-99,100)]
  
  set.seed(123)
  partition_index = createDataPartition(as.factor(gmfcs), p = .7, 
                                 list = FALSE, 
                                 times = 10)
  
  # fPCA-SVM
  
  print('----------- SVM modeling ----------')
  
  pcs_list = vector('list',6)
  names(pcs_list) = columns[1:6]
  
  for (col in columns[1:6]){
    pcs_list[[col]] = fpc.score(to_matrix(data,col),100,3)
    pcs_list[[col]]$gmfcs = as.factor(gmfcs)
  }
  
  svm_perf = perf_matrix()
  
  for (col in columns[1:6]){
    
    print(paste(col, 'started'))
    
    col_data = pcs_list[[col]]
    for (i in seq(1,10)){
      
      print(paste(i,'th cv done.', sep = ''))
      
      cm = SVM_CV(col_data, 5, manual = TRUE, manual_index = partition_index[,i])
      
      svm_perf[[col]][[i]] = cm
    }
  }
  
  
  ### FLR 
  
  print('---------- FLR modeling ------------')
  
  flr_perf = perf_matrix()
  
  coef_list = vector('list',6)
  names(coef_list) = columns[1:6]
  
  basis_list = vector('list',6)
  names(basis_list) = columns[1:6]
  
  for (col in columns[1:6]){
    coef_list[[col]] = vector('list',10)
    basis_list[[col]] = matrix(nrow = 10, ncol = 1)
    colnames(basis_list[[col]]) = c('n_basis')
  }
  
  for (col in columns[1:6]){
    col_data = t(to_matrix(data,col))
    
    print(paste(col, 'started'))
    
    for (i in seq(1,10)){
      
      
      result = FLR_CV(col_data, gmfcs, 5, manual = TRUE, manual_index = partition_index[,i])
      
      coef_list[[col]][[i]] = result$coef
      
      flr_perf[[col]][[i]] = result$cm
      
      basis_list[[col]][i,] = result$n_basis
      
      print(paste(i,'th cv done.', sep = ''))
    }
  }
  
  
  ### SFLDA
  
  print('---------- SFLDA modeling ------------')
  
  SFLDA_perf = perf_matrix()
  
  beta_list = vector('list',6)
  names(beta_list) = columns[1:6]
  
  param_list = vector('list',6)
  names(param_list) = columns[1:6]
  
  for (col in columns[1:6]){
    beta_list[[col]] = matrix(nrow = 100, ncol = 10)
    param_list[[col]] = matrix(nrow = 10, ncol = 2)
    colnames(param_list[[col]]) = c('tau','lambda')
  }
  
  for (col in columns[1:6]){
    
    print(paste(col, 'started'))
    
    col_data = t(to_matrix(data,col))
    df = data.frame(col_data)
    df['group'] = gmfcs
    
    for (i in seq(1:10)){
      result = SFLDA_CV(df, 5, manual = TRUE, manual_index = partition_index[,i])
      
      beta_list[[col]][,i] = result$b
      
      param_list[[col]][i,] = c(result$tau, result$lambda)
      
      cm = confusionMatrix(reference = as.factor(result$yte), data = as.factor(result$class), positive = '2')
      
      SFLDA_perf[[col]][[i]] = cm
      
      print(paste(i,'th cv done.', sep = ''))
    }
  }
  
  
  return(list('svm' = svm_perf, 'flr' = flr_perf,'coef_list' = coef_list ,'SFLDA' = SFLDA_perf, 
              'beta_list' = beta_list, 'param_list' = param_list, 'basis_list' = basis_list))
}

### normal vs patient
data_patient = data

data_patient[data_patient$gmfcs >= 1,'gmfcs'] = 1
rownames(data_patient) = NULL

comparison_patient = model_comparsion(data_patient)

save(comparison_patient, file = "./Paper/data/comparison_patient_3.RData")

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

comparison_01 = model_comparsion(data_01)
save(comparison_01, file = "./Paper/data/comparison_01_3.RData")

# 1 vs 2

data_12 = pick_two(data,1,2)
rownames(data_12) = NULL

comparison_12 = model_comparsion(data_12)
save(comparison_12, file = "./Paper/data/comparison_12_3.RData")

# 2 vs 3

data_23 = pick_two(data,2,3)
rownames(data_23) = NULL

comparison_23 = model_comparsion(data_23)
save(comparison_23, file = "./Paper/data/comparison_23_3.RData")

calculate_metric = function(lst){
  
  final_list = list()
  
  for (md in c('svm','flr','SFLDA')){
  
    metric = matrix(nrow = 6, ncol = 6)
    colnames(metric) = columns[1:6]
    rownames(metric) = c('acc','acc_std','fnr','fnr_std','for','for_std')
    
    for (col in columns[1:6]){
      acc = c()
      fnr = c()
      for_ = c()
      for (i in seq(1,10)){
        acc = c(acc, lst[[md]][[col]][[i]]$overall['Accuracy'])
        fnr = c(fnr, 1 - lst[[md]][[col]][[i]]$byClass['Sensitivity'])
        for_ = c(for_, 1 - lst[[md]][[col]][[i]]$byClass['Neg Pred Value'])
      }
      
      metric['acc',col] = mean(acc,na.rm = TRUE)
      metric['acc_std',col] = sd(acc, na.rm = TRUE)
      
      metric['fnr',col] = mean(fnr, na.rm = TRUE)
      metric['fnr_std',col] = sd(fnr, na.rm = TRUE)
      
      metric['for',col] = mean(for_, na.rm = TRUE)
      metric['for_std',col] = sd(for_, na.rm = TRUE)
      
      
    }
    
    final_list[[md]] = metric
      
  }
  
  return(final_list)
}

to_table = function(values, metric){
  tab = matrix(nrow = 18, ncol = 1)
  tab_std = matrix(nrow = 18, ncol = 1)
  
  metric_std = paste(metric,'std',sep = '_')
  
  for (i in seq(1,3)){
    tab[i,1] = values[[i]][metric,1]
    tab[i+3,1] = values[[i]][metric,3]
    tab[i+6,1] = values[[i]][metric,5]
    
    tab[i+9,1] = values[[i]][metric,2]
    tab[i+12,1] = values[[i]][metric,4]
    tab[i+15,1] = values[[i]][metric,6]
    
  }
  
  for (i in seq(1,3)){
    tab_std[i,1] = values[[i]][metric_std,1]
    tab_std[i+3,1] = values[[i]][metric_std,3]
    tab_std[i+6,1] = values[[i]][metric_std,5]
    
    tab_std[i+9,1] = values[[i]][metric_std,2]
    tab_std[i+12,1] = values[[i]][metric_std,4]
    tab_std[i+15,1] = values[[i]][metric_std,6]
    
  }
  
  return(list(tab = tab, tab_std = tab_std))
}

combine_table = function(lst1, lst2, lst3, lst4, metric, value){
 
  result = round(cbind(to_table(calculate_metric(lst1),metric = metric)[[value]], to_table(calculate_metric(lst2),metric = metric)[[value]],
              to_table(calculate_metric(lst3),metric = metric)[[value]], to_table(calculate_metric(lst4),metric = metric)[[value]]), 2)
  
  return(result)
   
}

acc_ = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'acc', value = 'tab')
acc_std = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'acc', value = 'tab_std')

fnr_ = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'fnr', value = 'tab')
fnr_std = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'fnr', value = 'tab_std')

for_ = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'for', value = 'tab')
for_std = combine_table(comparison_01, comparison_12, comparison_23, comparison_patient, metric = 'for', value = 'tab_std')

acc_table = paste(as.character(acc_), paste('(',as.character(acc_std), ')',sep = ''))
dim(acc_table) = c(18,4)

fnr_table = paste(as.character(fnr_), paste('(',as.character(fnr_std), ')',sep = ''))
dim(fnr_table) = c(18,4)

for_table = paste(as.character(for_), paste('(',as.character(for_std), ')',sep = ''))
dim(for_table) = c(18,4)

write.csv(data.frame(acc_table), './Paper/data/acc_table.csv')
write.csv(data.frame(fnr_table), './Paper/data/fnr_table.csv')
write.csv(data.frame(for_table), './Paper/data/for_table.csv')

### Multiple

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

model_comparison_multi = function(data_multi){
  
  multi_gmfcs = data_multi$gmfcs[seq(1,nrow(data_multi)-99,100)]
  
  set.seed(1)
  partition_index = createDataPartition(as.factor(multi_gmfcs), p = .7, 
                                        list = FALSE, 
                                        times = 10)
  
  # fPCA-SVM 
  
  pcs_list = vector('list',6)
  names(pcs_list) = columns[1:6]
  
  for (col in columns[1:6]){
    pcs_list[[col]] = fpc.score(to_matrix(data_multi,col),100,3)
    pcs_list[[col]]$gmfcs = as.factor(multi_gmfcs)
  }
  
  svm_perf = multi_perf_matrix()
  
  for (col in columns[1:6]){
    col_data = pcs_list[[col]]
    for (i in seq(1,10)){
      cm = SVM_CV(col_data, 5, manual = TRUE, manual_index = partition_index[,i], positive = '3')
      
      svm_perf[['acc']][i,col] = cm$overall[1]
      svm_perf[['1-sen']][[i,col]] = multi_sen(cm$table)
      svm_perf[['1-npv']][[i, col]] = multi_npv(cm$table)
      
    }
  }
  
  
  # FLR 
  
  flr_perf = multi_perf_matrix()
  
  
  for (col in columns[1:6]){
    col_data = t(to_matrix(data_multi,col))
    for (i in seq(1,10)){
      result = FLR_CV(col_data, multi_gmfcs, 5, manual = TRUE, manual_index = partition_index[,i], positive = '3')
      
      cm = result$cm
      
      flr_perf[['acc']][i,col] = cm$overall[1]
      flr_perf[['1-sen']][i,col] = multi_sen(cm$table)
      flr_perf[['1-npv']][i,col] = multi_npv(cm$table)
    }
  }
  
  # SFLDA_multi
  
  SFLDA_perf = multi_perf_matrix()
  
  beta_list = vector('list',6)
  names(beta_list) = columns[1:6]
  
  for (col in columns[1:6]){
    beta_list[[col]] = vector('list',3)
    for (j in seq(1,3)){
      beta_list[[col]][[j]] = matrix(nrow = 100, ncol = 10)
    }
  }
  
  for (col in columns[1:6]){
    col_data = t(to_matrix(data_multi,col))
    df = data.frame(col_data)
    df['group'] = multi_gmfcs
    
    for (i in seq(1:10)){
      result = SFLDA_multi_CV(df, 5, manual = TRUE, manual_index = partition_index[,i])
      print(paste(i,'th cv done.', sep = ''))
      
      for (j in seq(1:2)){
        beta_list[[col]][[j]][,i] = result$beta[,j]
      }
      
      cm = confusionMatrix(reference = as.factor(result$yte), data = as.factor(result$class), positive = '3')
      
      SFLDA_perf[['acc']][i,col] = cm$overall[1]
      SFLDA_perf[['1-sen']][i,col] = multi_sen(cm$table)
      SFLDA_perf[['1-npv']][[i,col]] = multi_npv(cm$table)
      
    }
    print(paste(col, 'done'))
  }
  
  return(list('svm' = svm_perf, 'flr' = flr_perf, 'SFLDA' = SFLDA_perf, 'b_list' = beta_list))
}

comparison_multi = model_comparison_multi(data_multi)

save(comparison_multi, file = "./Paper/data/comparison_multi_3.RData")

multi_to_table = function(lst, metric){
  tab = matrix(nrow = 18, ncol = 1)
  tab_std = matrix(nrow = 18, ncol = 1)
  
  models = c('svm','flr','SFLDA')
  
  for (i in seq(1,3)){
    tab[i,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[1]
    tab[i+3,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[3]
    tab[i+6,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[5]
    
    tab[i+9,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[2]
    tab[i+12,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[4]
    tab[i+15,1] = colMeans(lst[[models[i]]][[metric]], na.rm = TRUE)[6]
    
  }
  
  for (i in seq(1,3)){
    tab_std[i,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[1]
    tab_std[i+3,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[3]
    tab_std[i+6,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[5]
    
    tab_std[i+9,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[2]
    tab_std[i+12,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[4]
    tab_std[i+15,1] = apply(lst[[models[i]]][[metric]],2,sd, na.rm = TRUE)[6]
    
  }
  
  return(list(tab = tab, tab_std = tab_std))
  
}

multi_final = function(lst, metric){
 
  string = as.character(round(multi_to_table(lst,metric)[['tab']],2))
  std_string = paste('(',as.character(round(multi_to_table(lst,metric)[['tab_std']],2)),')',sep = '')
  

  final = paste(string, std_string)
  dim(final) = c(18,1)
  
  return(final)
}

acc_final = multi_final(comparison_multi, 'acc')
fnr_final = multi_final(comparison_multi, '1-sen')
for_final = multi_final(comparison_multi, '1-npv')

write.csv(data.frame(acc_final),'./Paper/data/acc_multi.csv')
write.csv(data.frame(fnr_final),'./Paper/data/fnr_multi.csv')
write.csv(data.frame(for_final),'./Paper/data/for_multi.csv')


# total

acc_save = cbind(acc_table, acc_final)
fnr_save = cbind(fnr_table, fnr_final)
for_save = cbind(for_table, for_final)

write.csv(data.frame(acc_save),'./Paper/data/acc_save.csv')
write.csv(data.frame(fnr_save),'./Paper/data/fnr_save.csv')
write.csv(data.frame(for_save),'./Paper/data/for_save.csv')

get_best_param = function(lst, col){
  acc = c()
  for (i in seq(1:10)){
    acc = c(acc,lst$SFLDA[[col]][[i]]$overall[1])
  }
  max_cv = which(acc == max(acc))
  
  return(matrix(lst$param_list[[col]][max_cv,],ncol = 2))
}

get_b_mat = function(df,lst){
  b_matrix = matrix(nrow = 100, ncol = 6)
  colnames(b_matrix) = columns[1:6]
  for (col in columns[1:6]){
    sdelta = SDelta(t(to_matrix(df[df$gmfcs == 0,],col)),t(to_matrix(df[df$gmfcs == 1,],col)))
    param = get_best_param(lst, col)[1,]
    b_matrix[,col] = SFLDA_get_b(sdelta$S,sdelta$mu0,sdelta$mu1, param[1], param[2])
  }
  
  return(b_matrix)
}

b_matrix_01 = get_b_mat(data_01,comparison_01)
b_matrix_12 = get_b_mat(data_12,comparison_12)
b_matrix_23 = get_b_mat(data_23,comparison_23)
b_matrix_patient = get_b_mat(data_patient,comparison_patient)

for (i in seq(1,6)){
  b_mat  = data.frame(cbind(b_matrix_01[,i],b_matrix_12[,i],b_matrix_23[,i],seq(1,100)))
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
  
  ggsave(filename =paste('./Paper/Plots/Methods/SFLDA_beta_',columns[i],'.png',sep = ''), plot = p,width = 10.47, height = 6.98)
}

for (i in seq(1,6)){
  b_mat  = data.frame(cbind(b_matrix_patient[,i], seq(1,100)))
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
  
  ggsave(filename =paste('./Paper/Plots/Methods/SFLDA_beta_patient_',columns[i],'.png',sep = ''), plot = g,width = 10.47, height = 6.98)
}


