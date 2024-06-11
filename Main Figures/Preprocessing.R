#### Garner, Kind et al. 2024 ####
## preprocessing ##

# This script performs a couple of preprocessing steps and saves the output
# to disk. Only required to run once if the output of this script is missing
# or whenever there are changes to this script. "Main.R" loads the output of
# this script for further analysis. 

# The saved output objects are:

# mi1_l_ar
# mi1_r_ar

# mi1_l_trans
# mi1_r_trans

# m6_l_pca 
# m6_r_pca 

# m6_l_x_form
# m6_r_x_form


# load required packages --------------------------------------------------

library(tidyverse)
library(natverse)
library(fafbseg)
library(RANN)

# load custom functions ---------------------------------------------------

source("R/functions.R")


# load codex data (v783) --------------------------------------------------


visual_neuron_types <- read_csv("data/v783_codex/visual_neuron_types.csv", 
                                col_types = cols(root_id = col_character())) |> 
  mutate(type = ifelse(type %in% c("MeMe_e10"), "MeMeDRA", type)) 


classification <- read_csv("data/v783_codex/classification.csv",
                           col_types = cols(root_id = col_character()))

connections <- read_csv("data/v783_codex/connections.csv", 
                        col_types = cols(pre_root_id = col_character(),
                                         post_root_id = col_character()))

labels <- read_csv("data/v783_codex/labels.csv", 
                   col_types = cols(root_id = col_character(), 
                                    supervoxel_id = col_character()))

synapse_coordinates <- read_csv("data/v783_codex/synapse_coordinates.csv", 
                                col_types = cols(pre_root_id = col_character(), 
                                                 post_root_id = col_character())) |> 
  fill(c("pre_root_id", "post_root_id"), .direction = "down")

garner_kind_labels <- read_csv("data/Garner_Kind_labels.csv", 
                               col_types = cols(hemisphere = col_skip(), 
                                                seg_id = col_character())) |> 
  rename(hemisphere = post_switch_hemisphere,
         root_id = seg_id,
         label = type) |> 
  mutate(label = ifelse(label %in% c("MeTu3c_dorsal", "MeTu3c_ventral"), "MeTu3c", label)) 
  


project_labels <- garner_kind_labels |>
  dplyr::select(root_id, label)

# Mi1 medulla maps (l & r) ------------------------------------------------

# Creates medulla columns based on Mi1 neurons for the right and left optic 
# lobe. 

fb_layer <- c(0.000, 0.100, 0.200, 0.300, 0.350, 0.430, 0.540, 0.660, 0.730, 0.910, 1.000)
m7c_layer <- c(0.000,  0.082, 0.262, 0.361, 0.459, 0.541, 0.623, 0.672, 0.762, 0.910, 1.000)
fafb_metu_paper_layer <- c(-0.03930698,  0.05537918,  0.17146181,  0.30822117,
                           0.34041553, 0.43168889,  0.50147245,  0.63131141, 
                           0.75365122,  0.92391202,  1.02194512)

# 0.000 0.089 0.199 0.327 0.358 0.444 0.510 0.632 0.747 0.908 1.000

names(fb_layer) <- paste("m", 0:10, sep = "")
names(m7c_layer) <- paste("m", 0:10, sep = "")
names(fafb_metu_paper_layer) <- paste("m", 0:10, sep = "")

layer_type <- fafb_metu_paper_layer #m7c_layer # or fb_layer

ly_ld <- c("m1" = 0.008036097,
          "m2" = 0.113420492,
          "m3" = 0.239841492,
          "m4" = 0.324318354,
          "m5" = 0.386052213,
          "m6" = 0.466580673,
          "m7" = 0.566391929,
          "m8" = 0.692481310,
          "m9" = 0.838781619,
          "m10" = 0.972928571)


mi1_l_id <- labels |> 
  left_join(classification, by = "root_id") |> 
  filter( side == "left" & label == "Mi1" & user_name %in% c("Emil Kind", "Lucy Houghton")) |> 
  rename(type = label) |> 
  select(root_id, type, side) |> 
  arrange(root_id)

mi1_l_syn <- synapse_coordinates |> 
  fill(c("pre_root_id", "post_root_id"), .direction = "down") |> 
  filter(pre_root_id %in% mi1_l_id$root_id | post_root_id %in% mi1_l_id$root_id)|> 
  mutate(prepost = case_when(pre_root_id == post_root_id ~ "auto",
                             pre_root_id %in% mi1_l_id$root_id ~ "pre",
                             post_root_id %in% mi1_l_id$root_id ~ "post",
                             TRUE ~ "Error"),
         root_id = case_when(prepost  == "pre" ~ pre_root_id,
                             prepost == "post" ~ post_root_id,
                             prepost == "auto" ~ pre_root_id,
                             TRUE ~ "Error")) |> 
  arrange(root_id)

mi1_l_map <- data.frame()
xyz <- data.frame()

for (i in 1:length(unique(mi1_l_syn$root_id))) {
  
  xyz <- mi1_l_syn |> 
    filter(root_id == unique(mi1_l_syn$root_id)[[i]]) |> 
    select(x,y,z)
  
  N <- nrow(xyz) 
  mean_xyz <- apply(xyz, 2, mean)
  xyz_pca   <- princomp(xyz) 
  dirVector <- xyz_pca$loadings[, 1]   # PC1
  xyz_fit <- matrix(rep(mean_xyz, each = N), ncol=3) + xyz_pca$score[, 1] %*% t(dirVector) 
  
  minx <- quantile(xyz_fit[,1], c(0.03))
  maxx <- quantile(xyz_fit[,1], c(0.97))
  
  minindex <- which(abs(xyz_fit[,1] - minx) == (min(abs(xyz_fit[,1] - minx), na.rm = TRUE)))
  maxindex <- which(abs(xyz_fit[,1] - maxx) == (min(abs(xyz_fit[,1] - maxx), na.rm = TRUE)))
  
  dat <- data.frame(as.list(xyz_fit[minindex[1],]),as.list(xyz_fit[maxindex[1],]))
  dat <- dat %>% rename(m0_x = x, m0_y = y, m0_z = z,
                        m10_x = x.1, m10_y = y.1, m10_z = z.1)
  
  dat$name <- unique(mi1_l_syn$root_id)[[i]]
  mi1_l_map <- rbind(mi1_l_map, dat)
  
}

topbottom_l <- mi1_l_map

bottom_l <- topbottom_l[,4:6] %>% rename("x" = "m10_x", "y" = "m10_y", "z" = "m10_z")
top_l <- topbottom_l[,1:3] %>% rename("x" = "m0_x", "y" = "m0_y", "z" = "m0_z")

mi1_l_map[,c('m1_x','m1_y','m1_z')] <- (layer_type['m1'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m2_x','m2_y','m2_z')] <- (layer_type['m2'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m3_x','m3_y','m3_z')] <- (layer_type['m3'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m4_x','m4_y','m4_z')] <- (layer_type['m4'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m5_x','m5_y','m5_z')] <- (layer_type['m5'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m6_x','m6_y','m6_z')] <- (layer_type['m6'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m7_x','m7_y','m7_z')] <- (layer_type['m7'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m8_x','m8_y','m8_z')] <- (layer_type['m8'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m9_x','m9_y','m9_z')] <- (layer_type['m9'] * (bottom_l - top_l) + (top_l))

mi1_l_map[,c('m0_x','m0_y','m0_z')] <- (layer_type['m0'] * (bottom_l - top_l) + (top_l))
mi1_l_map[,c('m10_x','m10_y','m10_z')] <- (layer_type['m10'] * (bottom_l - top_l) + (top_l))

temp <- mi1_l_map %>% 
  pivot_longer(
    cols = !name,
    names_to = c("layer", "axis"),
    names_sep = "_",
    values_to = "coordinate") %>% 
  dplyr::select(name, axis, layer, coordinate) %>% 
  arrange(layer,axis)

temp_names <- lapply(temp[,1:3], unique)

mi1_l_ar <- array(temp$coordinate, dim=c(nrow(mi1_l_map),3,11), dimnames = temp_names)

save(mi1_l_ar, file = "R/preprocessed_data/mi1_l_ar.RData")

# pch3d(mi1_l_ar[,,"m0"], radius=500,col="red",pch=16,alpha=1)
# pch3d(mi1_l_ar[,,"m10"], radius=500,col="blue",pch=16,alpha=1)
# for(i in 1:nrow(mi1_l_ar[,,"m10"])) segments3d(rbind(mi1_l_ar[i,,"m10"], mi1_l_ar[i,,"m0"]), col="green3")


mi1_r_id <- visual_neuron_types |> 
  filter(type %in% c("Mi1")) |> 
  select(root_id, type, side) |> 
  arrange(root_id)

mi1_r_syn <- synapse_coordinates |> 
  fill(c("pre_root_id", "post_root_id"), .direction = "down") |> 
  filter(pre_root_id %in% mi1_r_id$root_id | post_root_id %in% mi1_r_id$root_id) |> 
  mutate(prepost = case_when(pre_root_id == post_root_id ~ "auto",
                             pre_root_id %in% mi1_r_id$root_id ~ "pre",
                             post_root_id %in% mi1_r_id$root_id ~ "post",
                             TRUE ~ "Error"),
         root_id = case_when(prepost  == "pre" ~ pre_root_id,
                             prepost == "post" ~ post_root_id,
                             prepost == "auto" ~ pre_root_id,
                             TRUE ~ "Error")) |> 
  arrange(root_id)


mi1_r_map <- data.frame()
xyz <- data.frame()

for (i in 1:length(unique(mi1_r_syn$root_id))) {
  
  xyz <- mi1_r_syn |> 
    filter(root_id == unique(mi1_r_syn$root_id)[[i]]) |> 
    select(x,y,z)
  
  N <- nrow(xyz) 
  mean_xyz <- apply(xyz, 2, mean)
  xyz_pca   <- princomp(xyz) 
  dirVector <- xyz_pca$loadings[, 1]   # PC1
  xyz_fit <- matrix(rep(mean_xyz, each = N), ncol=3) + xyz_pca$score[, 1] %*% t(dirVector) 
  
  minx <- quantile(xyz_fit[,1], c(0.03))
  maxx <- quantile(xyz_fit[,1], c(0.97))
  
  minindex <- which(abs(xyz_fit[,1] - minx) == (min(abs(xyz_fit[,1] - minx), na.rm = TRUE)))
  maxindex <- which(abs(xyz_fit[,1] - maxx) == (min(abs(xyz_fit[,1] - maxx), na.rm = TRUE)))
  
  dat <- data.frame(as.list(xyz_fit[minindex[1],]),as.list(xyz_fit[maxindex[1],]))
  dat <- dat %>% rename(m10_x = x, m10_y = y, m10_z = z,
                        m0_x = x.1, m0_y = y.1, m0_z = z.1)
  
  dat$name <- unique(mi1_r_syn$root_id)[[i]]
  mi1_r_map <- rbind(mi1_r_map, dat)
  
}

topbottom_r <- mi1_r_map

bottom_r <- topbottom_r[,1:3] %>% rename("x" = "m10_x", "y" = "m10_y", "z" = "m10_z")
top_r <- topbottom_r[,4:6] %>% rename("x" = "m0_x", "y" = "m0_y", "z" = "m0_z")

mi1_r_map[,c('m1_x','m1_y','m1_z')] <- (layer_type['m1'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m2_x','m2_y','m2_z')] <- (layer_type['m2'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m3_x','m3_y','m3_z')] <- (layer_type['m3'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m4_x','m4_y','m4_z')] <- (layer_type['m4'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m5_x','m5_y','m5_z')] <- (layer_type['m5'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m6_x','m6_y','m6_z')] <- (layer_type['m6'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m7_x','m7_y','m7_z')] <- (layer_type['m7'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m8_x','m8_y','m8_z')] <- (layer_type['m8'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m9_x','m9_y','m9_z')] <- (layer_type['m9'] * (bottom_r - top_r) + (top_r))

mi1_r_map[,c('m0_x','m0_y','m0_z')] <- (layer_type['m0'] * (bottom_r - top_r) + (top_r))
mi1_r_map[,c('m10_x','m10_y','m10_z')] <- (layer_type['m10'] * (bottom_r - top_r) + (top_r))

temp <- mi1_r_map %>% 
  pivot_longer(
    cols = !name,
    names_to = c("layer", "axis"),
    names_sep = "_",
    values_to = "coordinate") %>% 
  dplyr::select(name, axis, layer, coordinate) %>% 
  arrange(layer,axis)

temp_names <- lapply(temp[,1:3], unique)

mi1_r_ar <- array(temp$coordinate, dim=c(nrow(mi1_r_map),3,11), dimnames = temp_names)

save(mi1_r_ar, file = "R/preprocessed_data/mi1_r_ar.RData")

# pch3d(mi1_r_ar[,,"m0"], radius=500,col="red",pch=16,alpha=1)
# pch3d(mi1_r_ar[,,"m10"], radius=500,col="blue",pch=16,alpha=1)
# for(i in 1:nrow(mi1_r_ar[,,"m10"])) segments3d(rbind(mi1_r_ar[i,,"m10"], mi1_r_ar[i,,"m0"]), col="green3")


# Transformation array ----------------------------------------------------



mi1_r_vec <- mi1_r_ar[,,"m5"] - mi1_r_ar[,,"m3"]
mi1_r_vec_norm <- vec_norm(mi1_r_vec)

m6_r_pca <- princomp(mi1_r_ar[,,"m6"]) 
save(m6_r_pca, file = "R/preprocessed_data/m6_r_pca.RData")

m6_r_x_form <- sweep(as.matrix(mi1_r_ar[,,"m6"]), 2, m6_r_pca$center) %*% m6_r_pca$loadings %>% as.data.frame() 
save(m6_r_x_form, file = "R/preprocessed_data/m6_r_x_form.RData")


mat <- matrix(nrow = 3, ncol = 3)
mi1_r_trans <- replicate(nrow(mi1_r_vec_norm), mat)

for (i in 1:nrow(mi1_r_vec_norm)) {
  mi1_r_trans[1,,i] <- mi1_r_vec_norm[i,]
  mi1_r_trans[2,,i] <- pracma::cross(mi1_r_vec_norm[i,], m6_r_pca$loadings[,2])
  mi1_r_trans[3,,i] <- - pracma::cross	(mi1_r_trans[2,,i],mi1_r_vec_norm[i,])
}

save(mi1_r_trans, file = "R/preprocessed_data/mi1_r_trans.RData")


# pch3d(mi1_r_ar[,,"m10"], radius=500,col="grey20",pch=16,alpha=1)
# for (i in 1:nrow(mi1_r_vec_norm)) {
#   segments3d(matrix(c(mi1_r_ar[i,,"m10"],(mi1_r_ar[i,,"m10"] - (mi1_r_trans[1,,i] * 5000))), byrow = T, nrow = 2), col="red", lwd=3)
#   segments3d(matrix(c(mi1_r_ar[i,,"m10"],(mi1_r_ar[i,,"m10"] - (mi1_r_trans[2,,i] * 5000))), byrow = T, nrow = 2), col="green", lwd=3)
#   segments3d(matrix(c(mi1_r_ar[i,,"m10"],(mi1_r_ar[i,,"m10"] - (mi1_r_trans[3,,i] * 5000))), byrow = T, nrow = 2), col="blue", lwd=3)
# }



mi1_l_vec <- mi1_l_ar[,,"m5"] - mi1_l_ar[,,"m3"]
mi1_l_vec_norm <- vec_norm(mi1_l_vec)

m6_l_pca <- princomp(mi1_l_ar[,,"m6"]) 
save(m6_l_pca, file = "R/preprocessed_data/m6_l_pca.RData")

m6_l_x_form <- sweep(as.matrix(mi1_l_ar[,,"m6"]), 2, m6_l_pca$center) %*% m6_l_pca$loadings %>% as.data.frame() 
save(m6_l_x_form, file = "R/preprocessed_data/m6_l_x_form.RData")

mat <- matrix(nrow = 3, ncol = 3)
mi1_l_trans <- replicate(nrow(mi1_l_vec_norm), mat)

for (i in 1:nrow(mi1_l_vec_norm)) {
  mi1_l_trans[1,,i] <- mi1_l_vec_norm[i,]
  mi1_l_trans[2,,i] <- pracma::cross(mi1_l_vec_norm[i,], m6_l_pca$loadings[,2])
  mi1_l_trans[3,,i] <- - pracma::cross	(mi1_l_trans[2,,i],mi1_l_vec_norm[i,])
}
save(mi1_l_trans, file = "R/preprocessed_data/mi1_l_trans.RData")



# pch3d(mi1_l_ar[,,"m10"], radius=500,col="grey20",pch=16,alpha=1)
# for (i in 1:nrow(mi1_l_vec_norm)) {
#   segments3d(matrix(c(mi1_l_ar[i,,"m10"],(mi1_l_ar[i,,"m10"] - (mi1_l_trans[1,,i] * 5000))), byrow = T, nrow = 2), col="red", lwd=3)
#   segments3d(matrix(c(mi1_l_ar[i,,"m10"],(mi1_l_ar[i,,"m10"] - (mi1_l_trans[2,,i] * 5000))), byrow = T, nrow = 2), col="green", lwd=3)
#   segments3d(matrix(c(mi1_l_ar[i,,"m10"],(mi1_l_ar[i,,"m10"] - (mi1_l_trans[3,,i] * 5000))), byrow = T, nrow = 2), col="blue", lwd=3)
# }



# Process large Synapse Table ---------------------------------------------

column_names <- c("id", "pre_pt_root_id", "post_pt_root_id", "connection_score", "cleft_score", 
                  "gaba", "ach", "glut", "oct", "ser", "da", "pre_pt_supervoxel_id", 
                  "post_pt_supervoxel_id", "neuropil", "post_pt_position_x", 
                  "post_pt_position_y", "post_pt_position_z", "pre_pt_position_x", 
                  "pre_pt_position_y", "pre_pt_position_z")

filter_function <- function(df) {
  df |> 
    filter(pre_pt_root_id %in% project_labels$root_id | post_pt_root_id %in% project_labels$root_id) 
}

process_chunk <- function(chunk, pos) {
  filtered_chunk <- filter_function(chunk)
  if (nrow(filtered_chunk) > 0) {
    write_csv(filtered_chunk, "data/project_synapse_table.csv", append = TRUE)
  }
}

read_csv_chunked("data/v783_codex/Cloud_SQL_fw_mat783_1_synapses_all_valid.csv", 
                 callback = DataFrameCallback$new(process_chunk), 
                 col_names = column_names, 
                 col_types = cols(id = col_double(), 
                                  pre_pt_root_id = col_character(), 
                                  post_pt_root_id = col_character(), 
                                  connection_score = col_double(),
                                  cleft_score = col_double(), 
                                  gaba = col_double(), 
                                  ach = col_double(), 
                                  glut = col_double(), 
                                  oct = col_double(), 
                                  ser = col_double(), 
                                  da = col_double(),
                                  pre_pt_supervoxel_id = col_character(), 
                                  post_pt_supervoxel_id = col_character(), 
                                  neuropil = col_character(),
                                  post_pt_position_x = col_double(), 
                                  post_pt_position_y = col_double(), 
                                  post_pt_position_z = col_double(), 
                                  pre_pt_position_x = col_double(), 
                                  pre_pt_position_y = col_double(), 
                                  pre_pt_position_z = col_double()),
                 chunk_size = 1000000)

project_synapse_table <- read_csv ("data/project_synapse_table.csv", 
                 col_names = column_names, 
                 col_types = cols(id = col_double(), 
                                  pre_pt_root_id = col_character(), 
                                  post_pt_root_id = col_character(), 
                                  connection_score = col_double(),
                                  cleft_score = col_double(), 
                                  gaba = col_double(), 
                                  ach = col_double(), 
                                  glut = col_double(), 
                                  oct = col_double(), 
                                  ser = col_double(), 
                                  da = col_double(),
                                  pre_pt_supervoxel_id = col_character(), 
                                  post_pt_supervoxel_id = col_character(), 
                                  neuropil = col_character(),
                                  post_pt_position_x = col_double(), 
                                  post_pt_position_y = col_double(), 
                                  post_pt_position_z = col_double(), 
                                  pre_pt_position_x = col_double(), 
                                  pre_pt_position_y = col_double(), 
                                  pre_pt_position_z = col_double()))





# subset project_synapse_table  for relevant neuron types -----------------

## (MeTu, ER, TuTu, TuBu, ExR1)

metu_r <- garner_kind_labels |> 
  filter(class == "MeTu" & hemisphere == "R") |> 
  arrange(root_id)

metu_r_input <- project_synapse_table |> 
  filter(post_pt_root_id %in% metu_r$root_id)|> 
  left_join(metu_r[,c("label","root_id")], by = c("post_pt_root_id" = "root_id")) |> 
  left_join(visual_neuron_types, by = c("pre_pt_root_id"="root_id")) |> 
  rename(metu_type = label,
         pre_type = type) |> 
  left_join(project_labels, by = c("pre_pt_root_id"="root_id"), suffix = c("", "_visual")) |> 
  mutate(pre_type = ifelse(is.na(pre_type) | pre_type %in% c("R7", "R8"), label, pre_type)) |> 
  select(-label)  
 
metu_r_output <- project_synapse_table |> 
  filter(pre_pt_root_id %in% metu_r$root_id)|> 
  left_join(metu_r[,c("label","root_id")], by = c("pre_pt_root_id" = "root_id")) |> 
  left_join(visual_neuron_types, by = c("post_pt_root_id"="root_id")) |> 
  rename(metu_type = label,
         post_type = type) |> 
  left_join(project_labels, by = c("post_pt_root_id"="root_id"), suffix = c("", "_visual")) |> 
  mutate(post_type = ifelse(is.na(post_type) | post_type %in% c("R7", "R8"), label, post_type)) |> 
  select(-label) 


metu_l <- garner_kind_labels |> 
  filter(class == "MeTu" & hemisphere == "L") |> 
  arrange(root_id)

metu_l_input <- project_synapse_table |> 
  filter(post_pt_root_id %in% metu_l$root_id)|> 
  left_join(metu_l[,c("label","root_id")], by = c("post_pt_root_id" = "root_id")) |> 
  left_join(visual_neuron_types, by = c("pre_pt_root_id"="root_id")) |> 
  rename(metu_type = label,
         pre_type = type) |> 
  left_join(project_labels, by = c("pre_pt_root_id"="root_id"), suffix = c("", "_visual")) |> 
  mutate(pre_type = ifelse(is.na(pre_type) | pre_type %in% c("R7", "R8"), label, pre_type)) |> 
  select(-label)  

metu_l_output <- project_synapse_table |> 
  filter(pre_pt_root_id %in% metu_l$root_id)|> 
  left_join(metu_l[,c("label","root_id")], by = c("pre_pt_root_id" = "root_id")) |> 
  left_join(visual_neuron_types, by = c("post_pt_root_id"="root_id")) |> 
  rename(metu_type = label,
         post_type = type) |> 
  left_join(project_labels, by = c("post_pt_root_id"="root_id"), suffix = c("", "_visual")) |> 
  mutate(post_type = ifelse(is.na(post_type) | post_type %in% c("R7", "R8"), label, post_type)) |> 
  select(-label) 


save(metu_r, file = "R/preprocessed_data/metu_r.RData")
save(metu_r_input, file = "R/preprocessed_data/metu_r_input.RData")
save(metu_r_output, file = "R/preprocessed_data/metu_r_output.RData")

save(metu_l, file = "R/preprocessed_data/metu_l.RData")
save(metu_l_input, file = "R/preprocessed_data/metu_l_input.RData")
save(metu_l_output, file = "R/preprocessed_data/metu_l_output.RData")


er_r_result <- neuron_syn_table("ER","R")
list2env(er_r_result, envir = .GlobalEnv)

er_l_result <- neuron_syn_table("ER","L")
list2env(er_l_result, envir = .GlobalEnv)

tutu_result <- neuron_syn_table("TuTu")
list2env(tutu_result, envir = .GlobalEnv)

tubu_r_result <- neuron_syn_table("TuBu","R")
list2env(tubu_r_result, envir = .GlobalEnv)

tubu_l_result <- neuron_syn_table("TuBu","L")
list2env(tubu_l_result, envir = .GlobalEnv)

exr1_l_result <- neuron_syn_table("ExR1","L", filter_column = "label")
list2env(exr1_l_result, envir = .GlobalEnv)

exr1_r_result <- neuron_syn_table("ExR1","R", filter_column = "label")
list2env(exr1_r_result, envir = .GlobalEnv)

save(er_r_id, file = "R/preprocessed_data/er_r_id.RData")
save(er_r_input, file = "R/preprocessed_data/er_r_input.RData")
save(er_r_output, file = "R/preprocessed_data/er_r_output.RData")

save(er_l_id, file = "R/preprocessed_data/er_l_id.RData")
save(er_l_input, file = "R/preprocessed_data/er_l_input.RData")
save(er_l_output, file = "R/preprocessed_data/er_l_output.RData")

save(tutu_id, file = "R/preprocessed_data/tutu_id.RData")
save(tutu_input, file = "R/preprocessed_data/tutu_input.RData")
save(tutu_output, file = "R/preprocessed_data/tutu_output.RData")

save(tubu_r_id, file = "R/preprocessed_data/tubu_r_id.RData")
save(tubu_r_input, file = "R/preprocessed_data/tubu_r_input.RData")
save(tubu_r_output, file = "R/preprocessed_data/tubu_r_output.RData")

save(tubu_l_id, file = "R/preprocessed_data/tubu_l_id.RData")
save(tubu_l_input, file = "R/preprocessed_data/tubu_l_input.RData")
save(tubu_l_output, file = "R/preprocessed_data/tubu_l_output.RData")

save(exr1_r_id, file = "R/preprocessed_data/exr1_r_id.RData")
save(exr1_r_input, file = "R/preprocessed_data/exr1_r_input.RData")
save(exr1_r_output, file = "R/preprocessed_data/exr1_r_output.RData")

save(exr1_l_id, file = "R/preprocessed_data/exr1_l_id.RData")
save(exr1_l_input, file = "R/preprocessed_data/exr1_l_input.RData")
save(exr1_l_output, file = "R/preprocessed_data/exr1_l_output.RData")