library(natverse)
library(tidyverse)
library(fafbseg)
library(googlesheets4)

mi1_test_set = read_sheet( ss="1PxDAEmIE61t0JN3PIq9Xch4oNWv4Z0NH27Ac5X5mYSQ", sheet = "Sheet1")
# mi1_test_set = mi1_test_set[1:111,]
# mi1_test_set = mi1_test_set %>% separate(XYZ, c("voxel_raw_x","voxel_raw_y","voxel_raw_z"), sep = ",",remove = F)
# 
# mi1_test_set$seg_id = flywire_xyz2id(matrix(c(as.numeric(mi1_test_set$voxel_raw_x),
#                                      as.numeric(mi1_test_set$voxel_raw_y),
#                                      as.numeric(mi1_test_set$voxel_raw_z)),
#                                    ncol = 3,
#                                    byrow = FALSE,
#                                    dimnames = list(NULL, c("X","Y","Z"))),
#                             rawcoords = TRUE,
#                             fast_root = TRUE)

# only run this if you one query all Mi1 synapses again. This can take ~1h!!
# otherwise just load the "Mi1_syn.RData" file.

# syn_data = flywire_ntpred(mi1_test_set$seg_id, cleft.threshold = 70) 
# Mi1_syn = syn_data

load("Mi1_syn.RData")

split_syn_data <- split(Mi1_syn,f=Mi1_syn$query)

Mi1_map = data.frame()
xyz = data.frame()

for (i in 1:length(split_syn_data)) {
  xyz = data.frame(x = split_syn_data[[names(split_syn_data[i])]][["pre_x"]],
                   y = split_syn_data[[names(split_syn_data[i])]][["pre_y"]],
                   z = split_syn_data[[names(split_syn_data[i])]][["pre_z"]])
  N <- nrow(xyz) 
  mean_xyz <- apply(xyz, 2, mean)
  xyz_pca   <- princomp(xyz) 
  dirVector <- xyz_pca$loadings[, 1]   # PC1
  xyz_fit <- matrix(rep(mean_xyz, each = N), ncol=3) + xyz_pca$score[, 1] %*% t(dirVector) 
  
  minx = quantile(xyz_fit[,1], c(0.01))
  maxx = quantile(xyz_fit[,1], c(0.99))
  
  minindex <- which(abs(xyz_fit[,1] - minx) == (min(abs(xyz_fit[,1] - minx), na.rm = TRUE)))
  maxindex <- which(abs(xyz_fit[,1] - maxx) == (min(abs(xyz_fit[,1] - maxx), na.rm = TRUE)))
  
  dat = data.frame(as.list(xyz_fit[minindex[1],]),as.list(xyz_fit[maxindex[1],]))
  dat = dat %>% rename(top_x = x, top_y = y, top_z = z,
                       bottom_x = x.1, bottom_y = y.1, bottom_z = z.1)
  dat$name = names(split_syn_data[i])
  Mi1_map = rbind(Mi1_map, dat)
}

Mi1_map_top = Mi1_map %>% select(top_x,top_y,top_z) %>% rename(x= top_x, y = top_y, z = top_z)
Mi1_map_bottom = Mi1_map %>% select(bottom_x,bottom_y,bottom_z) %>% rename(x= bottom_x, y = bottom_y, z = bottom_z)

nopen3d()
spheres3d(Mi1_map[,1:3], col = "red",radius = 1000)
spheres3d(Mi1_map[,4:6], col = "black",radius = 1000)
for(i in 1:length(Mi1_map$name)) segments3d(rbind(Mi1_map_top[i,], Mi1_map_bottom[i,]), col="green3")
rgl.viewpoint(fov=0,zoom=0.7, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))

rgl.snapshot("pngs/Mi1_map_m10.png")

htmlwidgets::saveWidget(rglwidget(), "Mi1_map.html")

w = rbind(Mi1_map_top, Mi1_map_bottom)
Mi1_msh = ashape3d(as.matrix(w), alpha = 60000) %>% as.mesh3d()
shade3d(Mi1_msh, alpha=0.1, col='gray')
wire3d(Mi1_msh)

Mi1_syn_sampled = sample_n(Mi1_syn, 50000)
Mi1_syn_sampled_msh = ashape3d(as.matrix(Mi1_syn_sampled[,3:5]), alpha = 60000) %>% as.mesh3d()
shade3d(Mi1_syn_sampled_msh, alpha=0.1, col='gray')
spheres3d(Mi1_syn_sampled[,3:5], col = "blue",radius = 250)
