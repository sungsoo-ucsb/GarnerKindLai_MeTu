library(natverse)
library(tidyverse)
library(fafbseg)
library(googlesheets4)
library(RColorBrewer)
library(factoextra)
library(nat.templatebrains)
library(nat.flybrains)
library(nat.jrcbrains)
library(cowplot)

# setup -------------------------------------------------------------------


# colors

#col_pal = brewer.pal(n = 6, name = "Set2")

gold =  "#ffd700" 
persian_green =  "#1b998b"
flickr_pink =  "#ed217c"
space_cadet =  "#2d3047"
light_salmon =  "#ff9b71"
pumpkin =  "#fa7921"

col_pal = c(gold,persian_green,flickr_pink,space_cadet,light_salmon,pumpkin)

MeTu_a_col = col_pal[1]
MeTu_pc_col = col_pal[2]
MeTu_pl_col = col_pal[3]

TuBu_R2abd_col = col_pal[1]
TuBu_R2c_col = col_pal[4]
TuBu_R3w_col = col_pal[5]
TuBu_R4d_col = col_pal[3]
TuBu_R4m_col = col_pal[2]
TuBu_R5_col = col_pal[6]

cluster_1_col = col_pal[1] 
cluster_2_col = col_pal[2]
cluster_3_col = col_pal[3]
cluster_4_col = col_pal[4]
cluster_5_col = col_pal[5]
cluster_6_col = col_pal[6]


# MeTu_a_col =  "#ffd700" #gold 
# MeTu_pc_col =  "#1b998b" # persian_green
# MeTu_pl_col =  "#ed217c" # flickr_pink


# reading date from google sheet ------------------------------------------


# TuBu and MeTu google sheet
n = read_sheet("1hVVnsaDqNxZhMtv_O24iWZ9aTMsceRR27TboLnLCXtA",sheet = "AVP")

# updating Seg_IDs
start_time <- Sys.time()
n = n %>% separate(XYZ, c("voxel_raw_x","voxel_raw_y","voxel_raw_z"), sep = ",",remove = F)
n$seg_id = flywire_xyz2id(matrix(c(as.numeric(n$voxel_raw_x),
                                   as.numeric(n$voxel_raw_y),
                                   as.numeric(n$voxel_raw_z)),
                                 ncol = 3,
                                 byrow = FALSE,
                                 dimnames = list(NULL, c("X","Y","Z"))),
                          rawcoords = TRUE,
                          fast_root = TRUE)
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

# check for duplications (should return empty dataframe)
n %>% group_by(seg_id) %>% filter(n()>1) %>%  arrange(seg_id)

n = n %>% 
  distinct(seg_id, .keep_all = TRUE) %>% 
  arrange(type_emil,type_dustin)

# write updated table back to google sheet
sheet_write(n, ss="1hVVnsaDqNxZhMtv_O24iWZ9aTMsceRR27TboLnLCXtA", sheet = "AVP")



# flywire URLs for MeTu & TuBu --------------------------------------------


# generate URL with three MeTu types (a,pl,pc) differently colored

metu_col = n %>% 
  filter(type_emil == "MeTu") %>% 
  mutate(col = case_when(type_dustin == "MeTu_a_l" | type_dustin == "MeTu_a_r" ~ MeTu_a_col,
                         type_dustin == "MeTu_pc_l" | type_dustin == "MeTu_pc_r" ~ MeTu_pc_col,
                         type_dustin == "MeTu_pl_l" | type_dustin == "MeTu_pl_r" ~ MeTu_pl_col)) %>% 
  select(seg_id,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

metu_link =ngl_add_colours(fw_url, metu_col)

# will open the link in browser 
browseURL(as.character(metu_link))

# generate URL with different tubu types (a la dustin) differently colored

col_tubu = n %>% 
  filter(type_emil == "TuBu") %>% 
  mutate(col = case_when(type_dustin == "TuBu_R2abd_l" | type_dustin == "TuBu_R2abd_r" ~ TuBu_R2abd_col,
                         type_dustin == "TuBu_R2c_l" | type_dustin == "TuBu_R2c_r" ~ TuBu_R2c_col,
                         type_dustin == "TuBu_R3w_l" | type_dustin == "TuBu_R3w_r" ~ TuBu_R3w_col,
                         type_dustin == "TuBu_R4d_l" | type_dustin == "TuBu_R4d_r" ~ TuBu_R4d_col,
                         type_dustin == "TuBu_R4m_l" | type_dustin == "TuBu_R4m_r" ~ TuBu_R4m_col,
                         type_dustin == "TuBu_R5_l" | type_dustin == "TuBu_R5_r" ~ TuBu_R5_col)) %>% 
  select(seg_id,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

tubu_link =ngl_add_colours(fw_url, col_tubu)

# will open the link in browser 
browseURL(as.character(tubu_link))


# data prep ---------------------------------------------------------------


metu_l = n %>% 
  filter(type_emil == "MeTu" & hemisphere == "L")

tubu_l = n %>% 
  filter(type_emil == "TuBu" & hemisphere == "L")

metu_r = n %>% 
  filter(type_emil == "MeTu" & hemisphere == "R")

tubu_r = n %>% 
  filter(type_emil == "TuBu" & hemisphere == "R")


metu_tubu_mat_l = flywire_adjacency_matrix(inputids = metu_l$seg_id, outputids = tubu_l$seg_id)
colnames(metu_tubu_mat_l) = tubu_l$type_dustin

metu_tubu_mat_comb_l = t(apply(t(metu_tubu_mat_l), 2, function(x) tapply(x, colnames(metu_tubu_mat_l), sum, na.rm = TRUE)))


metu_tubu_mat_r = flywire_adjacency_matrix(inputids = metu_r$seg_id, outputids = tubu_r$seg_id)
colnames(metu_tubu_mat_r) = tubu_l$type_dustin

metu_tubu_mat_comb_r = t(apply(t(metu_tubu_mat_r), 2, function(x) tapply(x, colnames(metu_tubu_mat_r), sum, na.rm = TRUE)))


# 6 clusters left
metu_l_6_cluster_info = hcut(metu_tubu_mat_comb_l, 
                             k = 6,
                             hc_func = c("hclust"), 
                             hc_method = "ward.D2", 
                             hc_metric = "euclidean",
                             stand = T)
fviz_dend(metu_l_6_cluster_info)

metu_l_6_cluster = enframe(metu_l_6_cluster_info$cluster)
table(metu_l_6_cluster$value)

metu_l_6_cluster_col = metu_l_6_cluster %>% 
  mutate(col = case_when(value == "1" ~ cluster_1_col,
                         value == "2" ~ cluster_2_col,
                         value == "3" ~ cluster_3_col,
                         value == "4" ~ cluster_4_col,
                         value == "5" ~ cluster_5_col,
                         value == "6" ~ cluster_6_col,
  )) %>% 
  select(name,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

metu_l_6_cluster_link =ngl_add_colours(fw_url, metu_l_6_cluster_col)
browseURL(as.character(metu_l_6_cluster_link))


# 6 clusters right
metu_r_6_cluster_info = hcut(metu_tubu_mat_comb_r, 
                             k = 6,
                             hc_func = c("hclust"), 
                             hc_method = "ward.D2", 
                             hc_metric = "euclidean",
                             stand = T)
fviz_dend(metu_r_6_cluster_info)

metu_r_6_cluster = enframe(metu_r_6_cluster_info$cluster)
table(metu_r_6_cluster$value)

metu_r_6_cluster_col = metu_r_6_cluster %>% 
  mutate(col = case_when(value == "1" ~ cluster_1_col,
                         value == "2" ~ cluster_2_col,
                         value == "3" ~ cluster_3_col,
                         value == "4" ~ cluster_4_col,
                         value == "5" ~ cluster_5_col,
                         value == "6" ~ cluster_6_col,
  )) %>% 
  select(name,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

metu_r_6_cluster_link =ngl_add_colours(fw_url, metu_r_6_cluster_col)
browseURL(as.character(metu_r_6_cluster_link))


# load mesh neurons -------------------------------------------------------


metu_l_6_cluster_n1 =  metu_l_6_cluster %>% 
  filter(value==1)
metu_l_6_cluster_n2 =  metu_l_6_cluster %>% 
  filter(value==2)
metu_l_6_cluster_n3 =  metu_l_6_cluster %>% 
  filter(value==3)
metu_l_6_cluster_n4 =  metu_l_6_cluster %>% 
  filter(value==4)
metu_l_6_cluster_n5 =  metu_l_6_cluster %>% 
  filter(value==5)
metu_l_6_cluster_n6 =  metu_l_6_cluster %>% 
  filter(value==6)

metu_l_6_cluster_n1_msh = read_cloudvolume_meshes(metu_l_6_cluster_n1$name)
metu_l_6_cluster_n2_msh = read_cloudvolume_meshes(metu_l_6_cluster_n2$name)
metu_l_6_cluster_n3_msh = read_cloudvolume_meshes(metu_l_6_cluster_n3$name)
metu_l_6_cluster_n4_msh = read_cloudvolume_meshes(metu_l_6_cluster_n4$name)
metu_l_6_cluster_n5_msh = read_cloudvolume_meshes(metu_l_6_cluster_n5$name)
metu_l_6_cluster_n6_msh = read_cloudvolume_meshes(metu_l_6_cluster_n6$name)

metu_a_r = n %>% 
  filter(type_dustin == "MeTu_a_r" )

metu_pc_r = n %>% 
  filter(type_dustin == "MeTu_pc_r" )

metu_pl_r = n %>% 
  filter(type_dustin == "MeTu_pl_r" )



metu_a_r_msh = read_cloudvolume_meshes(metu_a_r$seg_id)
metu_pc_r_msh = read_cloudvolume_meshes(metu_pc_r$seg_id)
metu_pl_r_msh = read_cloudvolume_meshes(metu_pl_r$seg_id)

TuBu_R2abd_l = n %>% 
  filter(type_dustin == "TuBu_R2abd_l") 

TuBu_R2c_l = n %>% 
  filter(type_dustin == "TuBu_R2c_l") 

TuBu_R3w_l = n %>% 
  filter(type_dustin == "TuBu_R3w_l") 

TuBu_R4d_l = n %>% 
  filter(type_dustin == "TuBu_R4d_l") 

TuBu_R4m_l = n %>% 
  filter(type_dustin == "TuBu_R4m_l" ) 

TuBu_R5_l = n %>% 
  filter(type_dustin == "TuBu_R5_l") 



TuBu_R2abd_r = n %>% 
  filter(type_dustin == "TuBu_R2abd_r") 

TuBu_R2c_r = n %>% 
  filter(type_dustin == "TuBu_R2c_r") 

TuBu_R3w_r = n %>% 
  filter(type_dustin == "TuBu_R3w_r") 

TuBu_R4d_r = n %>% 
  filter(type_dustin == "TuBu_R4d_r") 

TuBu_R4m_r = n %>% 
  filter(type_dustin == "TuBu_R4m_r" ) 

TuBu_R5_r = n %>% 
  filter(type_dustin == "TuBu_R5_r") 



TuBu_R2abd_l_msh = read_cloudvolume_meshes(TuBu_R2abd_l$seg_id)
TuBu_R2c_l_msh = read_cloudvolume_meshes(TuBu_R2c_l$seg_id)
TuBu_R3w_l_msh = read_cloudvolume_meshes(TuBu_R3w_l$seg_id)
TuBu_R4d_l_msh = read_cloudvolume_meshes(TuBu_R4d_l$seg_id)
TuBu_R4m_l_msh = read_cloudvolume_meshes(TuBu_R4m_l$seg_id)
TuBu_R5_l_msh = read_cloudvolume_meshes(TuBu_R5_l$seg_id)


TuBu_R2abd_r_msh = read_cloudvolume_meshes(TuBu_R2abd_r$seg_id)
TuBu_R2c_r_msh = read_cloudvolume_meshes(TuBu_R2c_r$seg_id)
TuBu_R3w_r_msh = read_cloudvolume_meshes(TuBu_R3w_r$seg_id)
TuBu_R4d_r_msh = read_cloudvolume_meshes(TuBu_R4d_r$seg_id)
TuBu_R4m_r_msh = read_cloudvolume_meshes(TuBu_R4m_r$seg_id)
TuBu_R5_r_msh = read_cloudvolume_meshes(TuBu_R5_r$seg_id)
# plot --------------------------------------------------------------------

m = matrix(c(0,0,0,
             1000000,0,0,
             0,0, 350000, 
             1000000, 0, 350000,
             0,400000,0,
             1000000,400000,0,
             0,400000, 350000, 
             1000000, 400000, 350000),ncol = 3,byrow = T)

set.seed(123)
nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))

plot3d(metu_a_r_msh, col = MeTu_a_col)
plot3d(metu_pc_r_msh, col = MeTu_pc_col)
plot3d(metu_pl_r_msh, col = MeTu_pl_col)

plot3d(metu_l_6_cluster_n1_msh, col = MeTu_a_col)
plot3d(metu_l_6_cluster_n2_msh, col = MeTu_a_col)
plot3d(metu_l_6_cluster_n3_msh, col = MeTu_a_col)
plot3d(metu_l_6_cluster_n4_msh, col = MeTu_pc_col)
plot3d(metu_l_6_cluster_n5_msh, col = MeTu_pc_col)
plot3d(metu_l_6_cluster_n6_msh, col = MeTu_pl_col)

rgl.snapshot("nov_meeting/all_MeTu_col_sub1_2_3.png")


nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))

plot3d(metu_a_r_msh, col = MeTu_a_col)


plot3d(metu_l_6_cluster_n1_msh, col = MeTu_a_col)
plot3d(metu_l_6_cluster_n2_msh, col = MeTu_a_col)
plot3d(metu_l_6_cluster_n3_msh, col = MeTu_a_col)
rgl.snapshot("nov_meeting/all_MeTu_col_sub1.png")

nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))


plot3d(metu_pc_r_msh, col = MeTu_pc_col)

plot3d(metu_l_6_cluster_n4_msh, col = MeTu_pc_col)
plot3d(metu_l_6_cluster_n5_msh, col = MeTu_pc_col)

rgl.snapshot("nov_meeting/all_MeTu_col_sub2.png")

nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))


plot3d(metu_pl_r_msh, col = MeTu_pl_col)

plot3d(metu_l_6_cluster_n6_msh, col = MeTu_pl_col)

rgl.snapshot("nov_meeting/all_MeTu_col_sub3.png")


nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))


plot3d(metu_l_6_cluster_n1_msh, col = col_pal[1])
plot3d(metu_l_6_cluster_n2_msh, col = col_pal[2])
plot3d(metu_l_6_cluster_n3_msh, col = col_pal[3])


rgl.snapshot("nov_meeting/MeTu_l_col_cluster_1_2_3.png")

nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))

plot3d(metu_l_6_cluster_n4_msh, col = col_pal[1])
plot3d(metu_l_6_cluster_n5_msh, col = col_pal[2])



rgl.snapshot("nov_meeting/MeTu_l_col_cluster_4_5.png")


nopen3d()
par3d('windowRect' = c(100,100,2000,1100))
shade3d(as.mesh3d(JFRC2NP.surf.fafb), alpha=0.05, col='gray', lit=T)
points3d(m, alpha=0, size=0) 
rgl.viewpoint(fov=0,zoom=0.5, userMatrix= rotationMatrix(0/180*pi,0,0,1) %*% rotationMatrix(180/180*pi,1,0,0))

plot3d(metu_l_6_cluster_n1_msh, col = col_pal[1])
plot3d(metu_l_6_cluster_n2_msh, col = col_pal[2])
plot3d(metu_l_6_cluster_n3_msh, col = col_pal[3])
plot3d(metu_l_6_cluster_n4_msh, col = col_pal[4])
plot3d(metu_l_6_cluster_n5_msh, col = col_pal[5])
plot3d(metu_l_6_cluster_n6_msh, col = col_pal[6])

rgl.snapshot("nov_meeting/all_l_MeTu_col_all_cluster.png")

# nt_prediction -----------------------------------------------------------

# cluster 1
neuron_nt = flywire_ntpred(metu_l_6_cluster_n1$name)

neuron_nt_sum = neuron_nt %>% 
  group_by(query) %>% 
  count(top.nt) %>% 
  replace(is.na(.), 0) %>% 
  mutate(percentage = round(n / sum(n)*100,1))

neuron_nt_means <- aggregate(percentage ~  top.nt, neuron_nt_sum, function(x) mean = round(mean(x),1))

ggplot(neuron_nt_sum, aes(y = percentage, x = top.nt, col = top.nt)) + 
  geom_line(aes(group = query), col = "grey", alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  ylim(0,100) +
  scale_x_discrete(labels=c("acetylcholine" = "ACh", 
                            "dopamine" = "DA",
                            "gaba" = "GABA", 
                            "glutamate" = "Glu", 
                            "octopamine" = "Oct",
                            "serotonin" = "5-HT")) +
  geom_text(data = neuron_nt_means, aes(label = percentage, y = 100)) +
  ggtitle("Cluster_1") +
  theme_cowplot(12) +
  theme(axis.title.x = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5)) 

ggsave("nov_meeting/cluster_1_nt_pred.png", width = 4, height = 4)

# cluster 2
neuron_nt = flywire_ntpred(metu_l_6_cluster_n2$name)

neuron_nt_sum = neuron_nt %>% 
  group_by(query) %>% 
  count(top.nt) %>% 
  replace(is.na(.), 0) %>% 
  mutate(percentage = round(n / sum(n)*100,1))

neuron_nt_means <- aggregate(percentage ~  top.nt, neuron_nt_sum, function(x) mean = round(mean(x),1))

ggplot(neuron_nt_sum, aes(y = percentage, x = top.nt, col = top.nt)) + 
  geom_line(aes(group = query), col = "grey", alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  ylim(0,100) +
  scale_x_discrete(labels=c("acetylcholine" = "ACh", 
                            "dopamine" = "DA",
                            "gaba" = "GABA", 
                            "glutamate" = "Glu", 
                            "octopamine" = "Oct",
                            "serotonin" = "5-HT")) +
  geom_text(data = neuron_nt_means, aes(label = percentage, y = 100)) +
  ggtitle("Cluster_2") +
  theme_cowplot(12) +
  theme(axis.title.x = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5)) 

ggsave("nov_meeting/cluster_2_nt_pred.png", width = 4, height = 4)

# cluster 3
neuron_nt = flywire_ntpred(metu_l_6_cluster_n3$name)

neuron_nt_sum = neuron_nt %>% 
  group_by(query) %>% 
  count(top.nt) %>% 
  replace(is.na(.), 0) %>% 
  mutate(percentage = round(n / sum(n)*100,1))

neuron_nt_means <- aggregate(percentage ~  top.nt, neuron_nt_sum, function(x) mean = round(mean(x),1))

ggplot(neuron_nt_sum, aes(y = percentage, x = top.nt, col = top.nt)) + 
  geom_line(aes(group = query), col = "grey", alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  ylim(0,100) +
  scale_x_discrete(labels=c("acetylcholine" = "ACh", 
                            "dopamine" = "DA",
                            "gaba" = "GABA", 
                            "glutamate" = "Glu", 
                            "octopamine" = "Oct",
                            "serotonin" = "5-HT")) +
  geom_text(data = neuron_nt_means, aes(label = percentage, y = 100)) +
  ggtitle("Cluster_3") +
  theme_cowplot(12) +
  theme(axis.title.x = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 0.5)) 

ggsave("nov_meeting/cluster_3_nt_pred.png", width = 4, height = 4)

# save as nrrd format for fiji --------------------------------------------


metu_l_6_cluster_n1_msh_jrc18f = xform_brain(metu_l_6_cluster_n1_msh, sample="FlyWire", reference = "JRC2018F")
metu_l_6_cluster_n2_msh_jrc18f = xform_brain(metu_l_6_cluster_n2_msh, sample="FlyWire", reference = "JRC2018F")
metu_l_6_cluster_n3_msh_jrc18f = xform_brain(metu_l_6_cluster_n3_msh, sample="FlyWire", reference = "JRC2018F")
metu_l_6_cluster_n4_msh_jrc18f = xform_brain(metu_l_6_cluster_n4_msh, sample="FlyWire", reference = "JRC2018F")
metu_l_6_cluster_n5_msh_jrc18f = xform_brain(metu_l_6_cluster_n5_msh, sample="FlyWire", reference = "JRC2018F")
metu_l_6_cluster_n6_msh_jrc18f = xform_brain(metu_l_6_cluster_n6_msh, sample="FlyWire", reference = "JRC2018F")

metu_a_r_msh_jrc18f = xform_brain(metu_a_r_msh, sample="FlyWire", reference = "JRC2018F")
metu_pc_r_msh_jrc18f = xform_brain(metu_pc_r_msh, sample="FlyWire", reference = "JRC2018F")
metu_pl_r_msh_jrc18f = xform_brain(metu_pl_r_msh, sample="FlyWire", reference = "JRC2018F")



metu_l_6_cluster_n1_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n1_msh_jrc18f), JRC2018F)
metu_l_6_cluster_n2_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n2_msh_jrc18f), JRC2018F)
metu_l_6_cluster_n3_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n3_msh_jrc18f), JRC2018F)
metu_l_6_cluster_n4_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n4_msh_jrc18f), JRC2018F)
metu_l_6_cluster_n5_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n5_msh_jrc18f), JRC2018F)
metu_l_6_cluster_n6_msh_jrc18f_im = as.im3d(xyzmatrix(metu_l_6_cluster_n6_msh_jrc18f), JRC2018F)

metu_a_r_msh_jrc18f_im = as.im3d(xyzmatrix(metu_a_r_msh_jrc18f), JRC2018F)
metu_pc_r_msh_jrc18f_im = as.im3d(xyzmatrix(metu_pc_r_msh_jrc18f), JRC2018F)
metu_pl_r_msh_jrc18f_im = as.im3d(xyzmatrix(metu_pl_r_msh_jrc18f), JRC2018F)

metu_a_l_msh_jrc18f_im = as.im3d(xyzmatrix(c(metu_l_6_cluster_n1_msh_jrc18f,
                                             metu_l_6_cluster_n2_msh_jrc18f,
                                             metu_l_6_cluster_n3_msh_jrc18f)), JRC2018F)

metu_pc_l_msh_jrc18f_im = as.im3d(xyzmatrix(c(metu_l_6_cluster_n4_msh_jrc18f,
                                              metu_l_6_cluster_n5_msh_jrc18f)), JRC2018F)

metu_pl_l_msh_jrc18f_im = as.im3d(xyzmatrix(c(metu_l_6_cluster_n6_msh_jrc18f)), JRC2018F)

write.im3d(metu_l_6_cluster_n1_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n1_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_l_6_cluster_n2_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n2_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_l_6_cluster_n3_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n3_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_l_6_cluster_n4_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n4_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_l_6_cluster_n5_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n5_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_l_6_cluster_n6_msh_jrc18f_im, "nov_meeting/metu_l_6_cluster_n6_msh_jrc18f.nrrd", dtype='byte')

write.im3d(metu_a_r_msh_jrc18f_im, "nov_meeting/metu_a_r_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_pc_r_msh_jrc18f_im, "nov_meeting/metu_pc_r_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_pl_r_msh_jrc18f_im, "nov_meeting/metu_pl_r_msh_jrc18f.nrrd", dtype='byte')


write.im3d(metu_a_l_msh_jrc18f_im, "nov_meeting/metu_a_l_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_pc_l_msh_jrc18f_im, "nov_meeting/metu_pc_l_msh_jrc18f.nrrd", dtype='byte')
write.im3d(metu_pl_l_msh_jrc18f_im, "nov_meeting/metu_pl_l_msh_jrc18f.nrrd", dtype='byte')


# connectivity -------------------------------------------------------------


tubu_r_in = flywire_partner_summary(tubu_r$seg_id, partners = "input")

tubu_r_in_sum = tubu_r_in %>% 
  group_by(pre_id) %>% 
  summarise(weight_sum = sum(weight))

tubu_r_in_no_metu = tubu_r_in_sum %>% 
  anti_join(metu_r,by = c("pre_id"="seg_id")) %>% 
  filter(pre_id != 0 & weight_sum > 5) %>% 
  arrange(desc(weight_sum))


tubu_r_in_no_metu_pre = flywire_ntpred(tubu_r_in_no_metu$pre_id)

new = tubu_r_in_no_metu_pre
new$post_id = as.character(new$post_id)
new = new %>% 
  inner_join(tubu_r,by = c("post_id"="seg_id")) %>% 
  filter(cleft_scores != 0 & pre_x < 98500*4) %>% 
  group_by(query) %>% 
  summarise(n = n())




tubu_l_in = flywire_partner_summary(tubu_l$seg_id, partners = "input")

tubu_l_in_sum = tubu_l_in %>% 
  group_by(pre_id) %>% 
  summarise(weight_sum = sum(weight))

tubu_l_in_no_metu = tubu_l_in_sum %>% 
  anti_join(metu_l,by = c("pre_id"="seg_id")) %>% 
  filter(pre_id != 0 & weight_sum > 2) %>% 
  arrange(desc(weight_sum))


tubu_l_in_no_metu_pre = flywire_ntpred(tubu_l_in_no_metu$pre_id)

new_l = tubu_l_in_no_metu_pre
new_l$post_id = as.character(new_l$post_id)
new_l = new_l %>% 
  inner_join(tubu_l,by = c("post_id"="seg_id")) %>% 
  filter(cleft_scores != 0 & pre_x > 165288*4) %>% 
  group_by(query) %>% 
  summarise(n = n())

# stuff -------------------------------------------------------------------



col_cl_n1 = col_cluster_6 %>% filter(col=="green")

cl_n1_link =ngl_add_colours(fw_url, col_cl_n1)
browseURL(as.character(cl_n1_link))


col_cl_n2 = col_cluster_6 %>% filter(col=="red")

cl_n2_link =ngl_add_colours(fw_url, col_cl_n2)
browseURL(as.character(cl_n2_link))


col_cl_n3 = col_cluster_6 %>% filter(col=="blue")

cl_n3_link =ngl_add_colours(fw_url, col_cl_n3)
browseURL(as.character(cl_n3_link))


col_cl_n4 = col_cluster_6 %>% filter(col=="yellow")

cl_n4_link =ngl_add_colours(fw_url, col_cl_n4)
browseURL(as.character(cl_n4_link))


col_cl_n5 = col_cluster_6 %>% filter(col=="magenta")

cl_n5_link =ngl_add_colours(fw_url, col_cl_n5)
browseURL(as.character(cl_n5_link))


col_cl_n6 = col_cluster_6 %>% filter(col=="orange")

cl_n6_link =ngl_add_colours(fw_url, col_cl_n6)
browseURL(as.character(cl_n6_link))



# MeTu_pl_l_url = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6469892333109248"
# MeTu_pl_r_url = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6055732562624512"
# MeTu_a_l_url  = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6609513935273984"
# MeTu_a_r_url  = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5116622024998912"
# MeTu_pc_l_url = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6265543459864576"
# MeTu_pc_r_url = "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/6322809701662720"
# 
# 
# MeTu_pl_l = data.frame(type = "MeTu_pl_l",seg_id = ngl_segments(MeTu_pl_l_url))
# MeTu_pl_r = data.frame(type = "MeTu_pl_r",seg_id = ngl_segments(MeTu_pl_r_url))
# MeTu_a_l = data.frame(type = "MeTu_a_l",seg_id = ngl_segments(MeTu_a_l_url))
# MeTu_a_r = data.frame(type = "MeTu_a_r",seg_id = ngl_segments(MeTu_a_r_url))
# MeTu_pc_l = data.frame(type = "MeTu_pc_l",seg_id = ngl_segments(MeTu_pc_l_url))
# MeTu_pc_r = data.frame(type = "MeTu_pc_r",seg_id = ngl_segments(MeTu_pc_r_url))
# 
# 
# metu = rbind(MeTu_pl_l, MeTu_pl_r, MeTu_a_l, MeTu_a_r, MeTu_pc_l, MeTu_pc_r)
# 
# metu$seg_id = flywire_latestid(metu$seg_id, method = "leaves")
# 
# n_joined = full_join(n,metu,by = "seg_id")

mat_l_sc = scale(mat_l)
mat_l_sc_dist = get_dist(mat_l_sc, method = "euclidean")
fviz_dist(mat_l_sc_dist)


# Compute hierarchical clustering
mat_hc <- hclust(mat_l_sc_dist, method = "ward.D2")

# Visualize
plot(mat_hc, cex = 0.5)
fviz_nbclust(mat_l_sc, kmeans, method = "gap_stat")


set.seed(123) # for reproducibility
km.res <- kmeans(mat_l_sc, 6, nstart = 25)
# Visualize
fviz_cluster(km.res, data = mat_l_sc, palette = "jco",
             ggtheme = theme_minimal())


res.hc <- hclust(dist(mat_l_sc),  method = "ward.D2")
fviz_dend(res.hc, cex = 0.5, k = 5, palette = "jco") 



s =hcut(mat_l, 
        k = 5,
        hc_func = c("hclust"), 
        hc_method = "ward.D2", 
        hc_metric = "euclidean",
        stand = T)
fviz_dend(s)

w = enframe(s$cluster)
table(w$value)

col_cluster = w %>% 
  mutate(col = case_when(value == "1" ~ "green",
                         value == "2" ~ "red",
                         value == "3" ~ "blue",
                         value == "4" ~ "yellow",
                         value == "5" ~ "pink",
  )) %>% 
  select(name,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

cluster_link =ngl_add_colours(fw_url, col_cluster)
browseURL(as.character(metu_link))


mat_l_cluster =hcut(mat_l_gp_c, 
                    k = 5,
                    hc_func = c("hclust"), 
                    hc_method = "ward.D2", 
                    hc_metric = "euclidean",
                    stand = T)
fviz_dend(mat_l_cluster)

clusters = enframe(mat_l_cluster$cluster)
table(clusters$value)

col_cluster = clusters %>% 
  mutate(col = case_when(value == "1" ~ "green",
                         value == "2" ~ "red",
                         value == "3" ~ "blue",
                         value == "4" ~ "yellow",
                         value == "5" ~ "magenta",
  )) %>% 
  select(name,col)

fw_url=with_segmentation('flywire', getOption('fafbseg.sampleurl'))

cluster_link =ngl_add_colours(fw_url, col_cluster)
browseURL(as.character(cluster_link))


exp_seg_id = "720575940618791558"

exp_input = flywire_partner_summary(exp_seg_id, partners = "input")

xyz_all = flywire_ntpred(exp_input$pre_id)

exp_input_xyz = xyz_all %>% 
  filter(post_id == exp_seg_id)

exp_input_xyz_fil = xyz_all %>% 
  filter(post_id == exp_seg_id & cleft_scores > 50)

exp_output = flywire_ntpred(exp_seg_id)

exp_output_xyz = exp_output %>% 
  filter(cleft_scores > 50)

exp_msh = read_cloudvolume_meshes(exp_seg_id)
