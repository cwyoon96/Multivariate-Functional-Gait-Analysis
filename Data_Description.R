### Data Description Section
setwd("/Volumes/GoogleDrive/My Drive/KAIST/HS Lab/Projects/CP_SFLDA")

library(progress)
library(ggplot2)
library(dplyr)
source('./Paper/R codes/common_func.R')

# load data
data = read.csv("./data/gait_full_210803.csv")

table(data$gmfcs) / 100

columns = colnames(data)
columns = columns[-(1:5)]

new_name = c('R HIP Flexibility', 'L HIP Flexibility', 'R KNEE Flexibility', 'L KNEE Flexibility',
             'R ANK Flexibility', 'L ANK Flexibility',' R Pelvis Tilt', 'L Plevis Tilt',
             'R HIP Rotation', 'L HIP Rotation', 'R Foot Orientation', 'L Foot Orientation')

names(new_name) = columns

# Delete GMFCS level 4
data = data[data$gmfcs != 4,]

# Draw mean of each gmfcs

gmfcs_0 = data[data$gmfcs == 0, ]
gmfcs_1 = data[data$gmfcs == 1, ]
gmfcs_2 = data[data$gmfcs == 2, ]
gmfcs_3 = data[data$gmfcs == 3, ]

mean_curve_plot = function(col){
  mean_0 = rowMeans(to_matrix(gmfcs_0,col))
  mean_1 = rowMeans(to_matrix(gmfcs_1,col))
  mean_2 = rowMeans(to_matrix(gmfcs_2,col))
  mean_3 = rowMeans(to_matrix(gmfcs_3,col))
  
  mean = data.frame(time = rep(seq(1:100),4), mean = c(mean_0,mean_1,mean_2,mean_3),
                    GMFCS = rep(c('level 0', 'level 1', 'level 2', 'level 3'),each = 100))
  g = ggplot(data = mean, aes(x = time, y = mean, group = GMFCS, color = GMFCS)) + geom_line(size = 3, aes(linetype = GMFCS)) +
    theme_bw() +
    theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size = 18), panel.background = element_rect(colour = "black", size=2)) +
    scale_x_continuous(limits=c(0, 101), expand = c(0, 0)) +
    labs(y = '', x= 'Gait cycle') +
    scale_color_manual(values = c("#4DAF4A","#377EB8","#984EA3","#E41A1C")) +
    scale_linetype_manual(values=c("solid","twodash","longdash","dotted"))
  
  return(g)
  
}

for (col in columns){
  g = mean_curve_plot(col)
  print(g)
  ggsave(filename = paste('./Paper/Plots/Data Description/mean_',col,'.png',sep = ''), 
         plot = g, width = 10.47, height = 6.98)
}



data = data[,c('subject','age','gmfcs','time','R_HIP_Flex_ANG','L_HIP_Flex_ANG','R_KNEE_Flex_ANG','L_KNEE_Flex_ANG',
            'R_ANK_Flex_ANG','L_ANK_Flex_ANG','R_Pelvis_Lat_Tilt', 'L_Pelvis_Lat_Tilt',
            'R_HIP_Rot_ANG','L_HIP_Rot_ANG','R_Foot_Orientation', 'L_Foot_Orientation')]

write.csv(data,"./Paper/data/preprocessed_data_2.csv",row.names = FALSE)
