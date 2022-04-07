to_matrix = function(df, col) {
  df$subject = as.character(df$subject)
  temp = df[,c('subject', col)]
  sub = temp[seq(1,nrow(temp) - 99,100),'subject']
  n_sub = length(sub)
  final_mat = matrix(nrow = 100, ncol = n_sub)
  
  colnames(final_mat) = sub
  
  for (s in sub) {
    temp2 = temp %>%
      filter(subject == s)
    
    final_mat[,s] = as.numeric(temp2[,col])
  }
  
  return(final_mat)
}