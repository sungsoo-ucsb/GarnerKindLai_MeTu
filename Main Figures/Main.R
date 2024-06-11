#### Garner, Kind et al. 2024 ####



msh_tolerance <- 5000

aotu_bb_x_min = 160500
aotu_bb_x_max = 174300
aotu_bb_y_min = 39000
aotu_bb_y_max = 47500
aotu_bb_z_min = 1300
aotu_bb_z_max = 2200

# load required packages --------------------------------------------------

library(tidyverse)
library(natverse)
library(fafbseg)
library(alphashape3d)
library(Morpho)
library(RANN)
library(cowplot)
library(rlang)
library(dendextend)
library(circlize)
library(ComplexHeatmap)
library(MASS)
library(png)

# load custom functions ---------------------------------------------------

source("R/functions.R")

# load codex data (v783) --------------------------------------------------


visual_neuron_types <- read_csv("data/v783_codex/visual_neuron_types.csv", 
                                col_types = cols(root_id = col_character()))|> 
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

# load preprocessed data --------------------------------------------------

folder_path <- "R/preprocessed_data"
rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
for (file in rdata_files) {
  load(file)
}

# load neuropil meshes ----------------------------------------------------
 
load("neuropil_meshes/JFRC2NP.surf.fafb.rda")

me_l_3_8_msh <- ashape3d(as.matrix(rbind(mi1_l_ar[,,"m3"], mi1_l_ar[,,"m8"])), alpha = 60000) %>% as.mesh3d()
me_l_0_10_msh <- ashape3d(as.matrix(rbind(mi1_l_ar[,,"m0"], mi1_l_ar[,,"m10"])), alpha = 60000) %>% as.mesh3d()

me_r_3_8_msh <- ashape3d(as.matrix(rbind(mi1_r_ar[,,"m3"], mi1_r_ar[,,"m8"])), alpha = 60000) %>% as.mesh3d()
me_r_0_10_msh <- ashape3d(as.matrix(rbind(mi1_r_ar[,,"m0"], mi1_r_ar[,,"m10"])), alpha = 60000) %>% as.mesh3d()


brain <- as.mesh3d(JFRC2NP.surf.fafb)
aotu_l_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="AOTU_R")
aotu_r_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="AOTU_L")
lobula_r_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LO_L")
lobula_l_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LO_R")
righthemi <- cutMeshPlane(brain,v1 = c(530000,0,0), v2 = c(530000,0,10), v3 = c(530000,10,0),keep.upper = F )
lefthemi <- cutMeshPlane(brain,v1 = c(530000,0,0), v2 = c(530000,0,10), v3 = c(530000,10,0),keep.upper = T )


# MeTu r ------------------------------------------------------------------

ii <- pointsinside(metu_r_input[,c("pre_pt_position_x","pre_pt_position_y","pre_pt_position_z")],me_r_3_8_msh, rval = 'distance') > - msh_tolerance
ii <- which(ii==TRUE)
metu_r_input_me <- metu_r_input[ii,]

## Get closest column for each synapse in medulla and get layer estimate 

v = seq(0.05,1, by = 0.05)
column_mat = approx3d(mi1_r_ar[,"x","m10"],
                      mi1_r_ar[,"y","m10"],
                      mi1_r_ar[,"z","m10"],
                      mi1_r_ar[,"x","m0"],
                      mi1_r_ar[,"y","m0"],
                      mi1_r_ar[,"z","m0"],
                      v = v)
rownames(column_mat) = rep(rownames(mi1_r_ar[,,"m10"]), length(v))

nn = nn2(column_mat, metu_r_input_me[,c("pre_pt_position_x","pre_pt_position_y","pre_pt_position_z")], k = 1)

metu_r_input_me  = cbind(metu_r_input_me,rownames(column_mat)[nn$nn.idx]) |>  
  rename("Mi1" = "rownames(column_mat)[nn$nn.idx]") |>  
  mutate(column =  match(Mi1, rownames(mi1_r_ar[,,"m10"])))

metu_r_me_input_neuron <- metu_r_input_me |> 
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id,pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  select(pre_pt_root_id, post_pt_root_id, metu_type, pre_type, syn_count)

metu_r_me_input_type  <- metu_r_me_input_neuron |> 
  group_by(post_pt_root_id, pre_type) |>  
  mutate(cell_count = n(), syn_count_sum = sum(syn_count)) |> 
  ungroup() |>  
  distinct(across(c(post_pt_root_id,pre_type)), .keep_all = T) |>  
  select(post_pt_root_id:pre_type, syn_count_sum, cell_count) |> 
  group_by(post_pt_root_id) |>  
  mutate(percentage = round(syn_count_sum/sum(syn_count_sum), 3)*100) |>  
  arrange(post_pt_root_id,desc(percentage)) 

metu_r_me_input_type_proofed <- metu_r_me_input_type |> 
  filter(post_pt_root_id %in% proofed)


metu_r_input_aotu <- metu_r_input |>  
  filter(post_pt_position_x >= aotu_bb_x_min *  4 & post_pt_position_x <= aotu_bb_x_max *  4 &
         post_pt_position_y >=  aotu_bb_y_min *  4 & post_pt_position_y <=  aotu_bb_y_max *  4 &
         post_pt_position_z >=   aotu_bb_z_min * 40 & post_pt_position_z <=   aotu_bb_z_max * 40)

metu_r_aotu_input_neuron = metu_r_input_aotu |>  
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id,pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  select(pre_pt_root_id, post_pt_root_id, metu_type, pre_type, syn_count)

metu_r_aotu_input_type  = metu_r_aotu_input_neuron |>  
  group_by(post_pt_root_id, pre_type) |>  
  mutate(cell_count = n(), syn_count_sum = sum(syn_count)) |> 
  ungroup() |>  
  distinct(across(c(post_pt_root_id,pre_type)), .keep_all = T) |>  
  select(post_pt_root_id:pre_type, syn_count_sum, cell_count) |> 
  group_by(post_pt_root_id) |>  
  mutate(percentage = round(syn_count_sum/sum(syn_count_sum), 3)*100) |>  
  arrange(post_pt_root_id,desc(percentage)) 

metu_r_aotu_input_type_proofed <- metu_r_aotu_input_type |> 
  filter(post_pt_root_id %in% proofed)



metu_r_output_aotu <- metu_r_output |>  
  filter(pre_pt_position_x >= aotu_bb_x_min *  4 & pre_pt_position_x <= aotu_bb_x_max *  4 &
           pre_pt_position_y >=  aotu_bb_y_min *  4 & pre_pt_position_y <=  aotu_bb_y_max *  4 &
           pre_pt_position_z >=   aotu_bb_z_min * 40 & pre_pt_position_z <=   aotu_bb_z_max * 40)

metu_r_aotu_output_neuron = metu_r_output_aotu |>  
  group_by(pre_pt_root_id, post_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(pre_pt_root_id, post_pt_root_id)), .keep_all = T) |> 
  arrange(pre_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("post_type"), ~replace(., is.na(.), "Unidentified")) |> 
  select(pre_pt_root_id, post_pt_root_id, metu_type, post_type, syn_count)

metu_r_aotu_output_type  = metu_r_aotu_output_neuron |>  
  group_by(pre_pt_root_id, post_type) |>  
  mutate(cell_count = n(), syn_count_sum = sum(syn_count)) |> 
  ungroup() |>  
  distinct(across(c(pre_pt_root_id,post_type)), .keep_all = T) |>  
  select(pre_pt_root_id, metu_type, post_type, syn_count_sum, cell_count) |> 
  group_by(pre_pt_root_id) |>  
  mutate(percentage = round(syn_count_sum/sum(syn_count_sum), 3)*100) |>  
  arrange(pre_pt_root_id,desc(percentage)) 

metu_r_aotu_output_type_proofed <- metu_r_aotu_output_type |> 
  filter(pre_pt_root_id %in% proofed)

# MeTu l ------------------------------------------------------------------


ii <- pointsinside(metu_l_input[,c("pre_pt_position_x","pre_pt_position_y","pre_pt_position_z")],me_l_3_8_msh, rval = 'distance') > - msh_tolerance
ii <- which(ii==TRUE)
metu_l_input_me <- metu_l_input[ii,]

## Get closest column for each synapse in medulla and get layer estimate 

v = seq(0.05,1, by = 0.05)
column_mat = approx3d(mi1_l_ar[,"x","m10"],
                      mi1_l_ar[,"y","m10"],
                      mi1_l_ar[,"z","m10"],
                      mi1_l_ar[,"x","m0"],
                      mi1_l_ar[,"y","m0"],
                      mi1_l_ar[,"z","m0"],
                      v = v)
rownames(column_mat) = rep(rownames(mi1_l_ar[,,"m10"]), length(v))

nn = nn2(column_mat, metu_l_input_me[,c("pre_pt_position_x","pre_pt_position_y","pre_pt_position_z")], k = 1)

metu_l_input_me  = cbind(metu_l_input_me,rownames(column_mat)[nn$nn.idx]) |>  
  rename("Mi1" = "rownames(column_mat)[nn$nn.idx]") |>  
  mutate(column =  match(Mi1, rownames(mi1_l_ar[,,"m10"])))

metu_l_me_input_neuron <- metu_l_input_me |> 
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id,pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  select(pre_pt_root_id, post_pt_root_id, metu_type, pre_type, syn_count)




### Top5 input for all MeTu_r

# ls <- vector("list", length = 10)
# ls2 <- vector("list", length = 10)
# names(ls) <- sort(unique(metu_r$label))
# k <- 1
# for (i in sort(unique(metu_r$label))) {
#   tmp <- metu_r_me_input_type |>  
#     filter(metu_type == i) |>  
#     select(post_pt_root_id, pre_type, percentage) %>% 
#     pivot_wider(names_from = pre_type, values_from = percentage, values_fill = 0)  |>  
#     pivot_longer(cols = !post_pt_root_id,names_to = "pre_type", values_to = "percentage") 
#   tmp2 <- tmp |> 
#     group_by(pre_type) |>  
#     summarise(mean_percentage = round(mean(percentage),1),
#               sd_percentage = round(sd(percentage),1),
#               min_percentage = min(percentage),
#               max_percentage = max(percentage)) |>  
#     arrange(desc(mean_percentage)) |>  
#     filter(pre_type != "Unidentified") |>  
#     slice(1:5) |>  
#     pull(pre_type)
#   ls[[k]] <- c(tmp2, "Unidentified")
#   ls2[[k]] <- tmp |> filter(pre_type %in% ls[[k]])
#   k <- k +1
# }

### Top5 input for all proofread MeTu_r neurons

ls_proofed <- vector("list", length = 10)
ls2_proofed <- vector("list", length = 10)
names(ls_proofed) <- sort(unique(metu_r$label))
k <- 1
for (i in sort(unique(metu_r$label))) {
  tmp <- metu_r_me_input_type_proofed |>  
    filter(metu_type == i) |>  
    select(post_pt_root_id, pre_type, percentage) %>% 
    pivot_wider(names_from = pre_type, values_from = percentage, values_fill = 0)  |>  
    pivot_longer(cols = !post_pt_root_id,names_to = "pre_type", values_to = "percentage") 
  tmp2 <- tmp |> 
    group_by(pre_type) |>  
    summarise(mean_percentage = round(mean(percentage),1),
              sd_percentage = round(sd(percentage),1),
              min_percentage = min(percentage),
              max_percentage = max(percentage)) |>  
    arrange(desc(mean_percentage)) |>  
    filter(pre_type != "Unidentified") |>  
    slice(1:5) |>  
    pull(pre_type)
  ls_proofed[[k]] <- c(tmp2, "Unidentified")
  ls2_proofed[[k]] <- tmp |> filter(pre_type %in% ls_proofed[[k]])
  k <- k +1
}




ls2 <- vector("list", length = 10)
#names(ls) <- sort(unique(metu_r$label))
k <- 1
for (i in sort(unique(metu_r$label))) {
  tmp <- metu_r_me_input_type |>  
    filter(metu_type == i) |>  
    select(post_pt_root_id, pre_type, percentage) %>% 
    pivot_wider(names_from = pre_type, values_from = percentage, values_fill = 0)  |>  
    pivot_longer(cols = !post_pt_root_id,names_to = "pre_type", values_to = "percentage") 
  ls2[[k]] <- tmp |> filter(pre_type %in% ls_proofed[[k]])
  k <- k +1
}
# Top5 medulla input ------------------------------------------------------

means <- ls2[[1]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl1 = ggplot(ls2[[1]], aes(y = percentage , x = factor(pre_type, level = ls[[1]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[1]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())


means <- ls2[[2]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl2a = ggplot(ls2[[2]], aes(y = percentage , x = factor(pre_type, level = ls[[2]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[2]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[3]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl2b = ggplot(ls2[[3]], aes(y = percentage , x = factor(pre_type, level = ls[[3]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[3]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[4]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl3a = ggplot(ls2[[4]], aes(y = percentage , x = factor(pre_type, level = ls[[4]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[4]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[5]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl3b = ggplot(ls2[[5]], aes(y = percentage , x = factor(pre_type, level = ls[[5]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[5]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[6]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl3c = ggplot(ls2[[6]], aes(y = percentage , x = factor(pre_type, level = ls[[6]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[6]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[7]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl4a = ggplot(ls2[[7]], aes(y = percentage , x = factor(pre_type, level = ls[[7]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[7]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[8]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl4b = ggplot(ls2[[8]], aes(y = percentage , x = factor(pre_type, level = ls[[8]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[8]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())

means <- ls2[[9]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl4c = ggplot(ls2[[9]], aes(y = percentage , x = factor(pre_type, level = ls[[9]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[9]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())


means <- ls2[[10]] %>%
  group_by(pre_type) %>%
  summarise(mean_value = mean(percentage),
            sd_value = sd(percentage))

pl4d = ggplot(ls2[[10]], aes(y = percentage , x = factor(pre_type, level = ls[[10]]))) +
  geom_jitter(width = 0.15, alpha=1, size = 0.5, stroke = 0, shape = 16) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", color = "red", width = 0.3, linewidth = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.45, linewidth = 0.1) +
  geom_text(data = means, aes(x = factor(pre_type, levels = ls[[10]]), y = 100, label = paste0(round(mean_value, 1), " %")), hjust = 1,
            size = 1.5, color = "black", angle = 90) +
  ylab("(%) medulla input") +
  ylim(-2,102)+
  # annotate("text", x=6, y=50, label= paste("n =",length(unique(metu1_top5_syn$query)))) + 
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
  cowplot::theme_cowplot(font_size = 6,
                         line_size = 0.3,
                         font_family = "") +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())


plot_grid(pl1, ncol = 1, align = "h")
ggsave("output/MeTu1_top5.pdf", units = "mm", width = 23, height = 36)
plot_grid(pl2a,pl2b, ncol = 2, align = "h")
ggsave("output/MeTu2_top5.pdf", units = "mm", width = 23*2, height = 36)
plot_grid(pl3a,pl3b,pl3c, ncol = 3, align = "h")
ggsave("output/MeTu3_top5.pdf", units = "mm", width = 23*3, height = 36)
plot_grid(pl4a,pl4b,pl4c,pl4d, ncol = 4, align = "h")
ggsave("output/MeTu4_top5.pdf", units = "mm", width = 23*4, height = 36)










pl_ls <- list()
for (i in 1:length(ls2)) {
  ls2_combined <- bind_rows(
    ls2_proofed[[i]]  |>  mutate(dataset = "proofed"),
    ls2[[i]] |> mutate(dataset = "all")
  )
  
  ls2_combined$dataset <- factor(ls2_combined$dataset, levels = c("proofed", "all"))
  
  means_combined <- ls2_combined |>
    group_by(pre_type, dataset) |>
    summarise(mean_value = mean(percentage), .groups = 'drop')
  
  overall_means <- means_combined |>
    group_by(pre_type) |>
    summarise(overall_mean = mean(mean_value)) |>
    arrange(desc(overall_mean))
  
  overall_means <- overall_means |>
    filter(pre_type != "Unidentified") |>
    bind_rows(tibble(pre_type = "Unidentified", overall_mean = NA))
  
  ls2_combined$pre_type <- factor(ls2_combined$pre_type, levels = overall_means$pre_type)
  means_combined$pre_type <- factor(means_combined$pre_type, levels = overall_means$pre_type)
  
  position <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)
  
  custom_colors <- c("proofed" = "black", "all" = "dimgrey")
  
  pl_ls[[i]] <- ggplot(ls2_combined, aes(y = percentage, x = pre_type, color = dataset)) +
    geom_jitter(position = position, alpha = 1, size = 0.5, stroke = 0, shape = 16, show.legend = FALSE) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", position = position_dodge(width = 0.8), aes(group = dataset), color = "red2" , width = 0.3, linewidth = 0.15, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "crossbar", position = position_dodge(width = 0.8), aes(group = dataset), color = "red2", width = 0.45, linewidth = 0.1, show.legend = FALSE) +
    geom_text(data = means_combined, aes(x = pre_type, color = dataset, y = 102, label = paste0(round(mean_value, 1), "%"), group = dataset),
              position = position_dodge(width = 0.8), hjust = 1, size = 1.5, angle = 90, show.legend = FALSE) +
    scale_color_manual(values = custom_colors) +
    ylab("(%) medulla input") +
    ylim(-2, 102) +
    scale_x_discrete(guide = guide_axis(angle = 90),
                     labels = function(x) ifelse(x == "MeMe; I only identified it, have not actively checked proofreading", "MeMe", x)) +
    cowplot::theme_cowplot(font_size = 6,
                           line_size = 0.3,
                           font_family = "") +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
}

plot_grid(pl_ls[[1]], ncol = 1, align = "h")
ggsave("output/MeTu1_top5_both.pdf", units = "mm", width = 23*2, height = 36*1.2)
plot_grid(pl_ls[[2]],pl_ls[[3]], ncol = 2, align = "h")
ggsave("output/MeTu2_top5_both.pdf", units = "mm", width = 23*2*2, height = 36*1.2)
plot_grid(pl_ls[[4]],pl_ls[[5]],pl_ls[[6]], ncol = 3, align = "h")
ggsave("output/MeTu3_top5_both.pdf", units = "mm", width = 23*3*2, height = 36*1.2)
plot_grid(pl_ls[[7]],pl_ls[[8]],pl_ls[[9]],pl_ls[[10]], ncol = 4, align = "h")
ggsave("output/MeTu4_top5_both.pdf", units = "mm", width = 23*4*2, height = 36*1.2)

plot_grid(pl_ls[[1]],pl_ls[[2]],pl_ls[[3]],
          pl_ls[[4]],pl_ls[[5]],pl_ls[[6]],
          pl_ls[[7]],pl_ls[[8]],pl_ls[[9]],
          pl_ls[[10]], ncol = 10, align = "h")

ggsave("output/MeTu_top5_both.pdf", units = "mm", width = 23*10*1.5, height = 36)




conn_matrix_plot(metutype = c("MeTu1"), 
                 posttype = c("TuBu08"), 
                 cluster_n = 1, 
                 file_name = "MeTu1_dend_matrix_proofed")
conn_matrix_plot(metutype = c("MeTu2a","MeTu2b"), 
                 posttype = c("TuBu01", "TuBu06"), 
                 cluster_n = 2,
                 file_name = "MeTu2_dend_matrix_proofed")
conn_matrix_plot(metutype = c("MeTu3a","MeTu3b","MeTu3c"), 
                 posttype = c("TuBu07", "TuBu09", "TuBu10"), 
                 cluster_n = 3,
                 file_name = "MeTu3_dend_matrix_proofed")
conn_matrix_plot(metutype = c("MeTu4a","MeTu4b","MeTu4c","MeTu4d"),
                 posttype = c("TuBu02", "TuBu03", "TuBu04", "TuBu05"),
                 cluster_n = 4,
                 file_name = "MeTu4_dend_matrix_proofed_4_cluster")
conn_matrix_plot(metutype = c("MeTu4a","MeTu4b","MeTu4c","MeTu4d"),
                 posttype = c("TuBu02", "TuBu03", "TuBu04", "TuBu05"),
                 cluster_n = 5,
                 file_name = "MeTu4_dend_matrix_proofed_5_cluster")
conn_matrix_plot(metutype = c("MeTu1","MeTu2a","MeTu2b",
                              "MeTu3a","MeTu3b","MeTu3c",
                              "MeTu4a","MeTu4b","MeTu4c",
                              "MeTu4d"), 
                 posttype = c("TuBu01","TuBu02","TuBu03",
                              "TuBu04", "TuBu05","TuBu06",
                              "TuBu07", "TuBu08", "TuBu09",
                              "TuBu10"),
                 cluster_n = 10,
                 file_name = "MeTu_dend_matrix_proofed_10_cluster")
conn_matrix_plot(metutype = c("MeTu1","MeTu2a","MeTu2b",
                              "MeTu3a","MeTu3b","MeTu3c",
                              "MeTu4a","MeTu4b","MeTu4c",
                              "MeTu4d"), 
                 posttype = c("TuBu01","TuBu02","TuBu03",
                              "TuBu04", "TuBu05","TuBu06",
                              "TuBu07", "TuBu08", "TuBu09",
                              "TuBu10"),
                 cluster_n = 11,
                 file_name = "MeTu_dend_matrix_proofed_11_cluster")



conn_matrix_plot(metutype = c("MeTu1"), 
                 posttype = c("TuBu08"), 
                 cluster_n = 1, 
                 population = T,
                 file_name = "MeTu1_dend_matrix")
conn_matrix_plot(metutype = c("MeTu2a","MeTu2b"), 
                 posttype = c("TuBu01", "TuBu06"), 
                 cluster_n = 2,
                 population = T,
                 file_name = "MeTu2_dend_matrix")
conn_matrix_plot(metutype = c("MeTu3a","MeTu3b","MeTu3c"), 
                 posttype = c("TuBu07", "TuBu09", "TuBu10"), 
                 cluster_n = 3,
                 population = T,
                 file_name = "MeTu3_dend_matrix")
conn_matrix_plot(metutype = c("MeTu4a","MeTu4b","MeTu4c","MeTu4d"),
                 posttype = c("TuBu02", "TuBu03", "TuBu04", "TuBu05"),
                 cluster_n = 4,
                 population = T,
                 file_name = "MeTu4_dend_matrix_4_cluster")
conn_matrix_plot(metutype = c("MeTu4a","MeTu4b","MeTu4c","MeTu4d"),
                 posttype = c("TuBu02", "TuBu03", "TuBu04", "TuBu05"),
                 cluster_n = 5,
                 population = T,
                 file_name = "MeTu4_dend_matrix_5_cluster")
conn_matrix_plot(metutype = c("MeTu4a","MeTu4b","MeTu4c","MeTu4d"),
                 posttype = c("TuBu02", "TuBu03", "TuBu04", "TuBu05"),
                 cluster_n = 6,
                 population = T,
                 file_name = "MeTu4_dend_matrix_6_cluster")
conn_matrix_plot(metutype = c("MeTu1","MeTu2a","MeTu2b",
                              "MeTu3a","MeTu3b","MeTu3c",
                              "MeTu4a","MeTu4b","MeTu4c",
                              "MeTu4d"), 
                 posttype = c("TuBu01","TuBu02","TuBu03",
                              "TuBu04", "TuBu05","TuBu06",
                              "TuBu07", "TuBu08", "TuBu09",
                              "TuBu10"),
                 cluster_n = 10,
                 population = T,
                 file_name = "MeTu_dend_matrix_10_cluster")
conn_matrix_plot(metutype = c("MeTu1","MeTu2a","MeTu2b",
                              "MeTu3a","MeTu3b","MeTu3c",
                              "MeTu4a","MeTu4b","MeTu4c",
                              "MeTu4d"), 
                 posttype = c("TuBu01","TuBu02","TuBu03",
                              "TuBu04", "TuBu05","TuBu06",
                              "TuBu07", "TuBu08", "TuBu09",
                              "TuBu10"),
                 cluster_n = 11,
                 population = T,
                 file_name = "MeTu_dend_matrix_11_cluster")







# df_mat_r = read_sheet("1hU_Wbb-uPclLBuE5s-BKmYN5snu7Io4jYKcfMdXgO2I", na = "NA",sheet = "left_lobe_v783")
# df_mat_l = read_sheet("1hU_Wbb-uPclLBuE5s-BKmYN5snu7Io4jYKcfMdXgO2I", na = "NA",sheet = "right_lobe_v783")

df_mat_r <- read_csv("data/df_mat_r.csv")
df_mat_l <- read_csv("data/df_mat_l.csv")

hex_grid_df_l = col_hex_grid(df_mat_l)
hex_grid_df_r = col_hex_grid(df_mat_r)

lens_Mi1_left <- read_csv("data/lens_Mi1_left.csv")
thetaphi_left <- read_csv("data/thetaphi_left.csv")
thetaphi_left$ind <- 1:nrow(thetaphi_left)

lens_Mi1_right <- read_csv("data/lens_Mi1_right.csv")
thetaphi_right <- read_csv("data/thetaphi_right.csv")
thetaphi_right$ind <- 1:nrow(thetaphi_right)



tmp = lens_Mi1_right |> 
  left_join(thetaphi_right, by = c("V1"="ind"))


eye_map_r = hex_grid_df_r |> 
  left_join(tmp, by = c("col" = "V2")) |> 
  select(col,theta,phi) |> 
  mutate(hemisphere = "r",
         all_id = paste(col,hemisphere, sep = "_"),
         theta = -theta+90)


tmp = lens_Mi1_left |> 
  left_join(thetaphi_left, by = c("V1"="ind"))

eye_map_l = hex_grid_df_l |> 
  left_join(tmp, by = c("col" = "V2")) |> 
  select(col,theta,phi) |> 
  mutate(hemisphere = "l",
         all_id = paste(col,hemisphere, sep = "_"),
         theta = -theta+90)

eye_map = rbind(eye_map_l,eye_map_r)

#  r ------------------------------------------------------------------


# adding exr1 to the er_r_input

exr1_r_input <- exr1_r_input |> 
  rename(er_type = exr1_type)

er_r_input_neuron <- er_r_input |> 
  bind_rows(exr1_r_input) |> 
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id, pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  dplyr::select(pre_pt_root_id, post_pt_root_id, er_type, pre_type, er_type, syn_count)


## TuBu

tubu_r_input_neuron <- tubu_r_input |> 
  mutate(pre_type = ifelse(pre_type %in% c("MeTu3c_dorsal", "MeTu3c_ventral"), "MeTu3c", pre_type)) |> 
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id, pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  dplyr::select(pre_pt_root_id, post_pt_root_id, tubu_type, pre_type, syn_count)


## TuTu 

tutu_input_neuron <- tutu_input |> 
  mutate(pre_type = ifelse(pre_type %in% c("MeTu3c_dorsal", "MeTu3c_ventral"), "MeTu3c", pre_type)) |> 
  group_by(post_pt_root_id, pre_pt_root_id) |>  
  mutate(syn_count = n()) |> 
  ungroup() |> 
  distinct(across(c(post_pt_root_id, pre_pt_root_id)), .keep_all = T) |> 
  arrange(post_pt_root_id, desc(syn_count)) |> 
  mutate_at(c("pre_type"), ~replace(., is.na(.), "Unidentified")) |> 
  dplyr::select(pre_pt_root_id, post_pt_root_id, tutu_type, pre_type, syn_count)


####

er_df <- er_r_input_neuron |> 
  group_by(post_pt_root_id) |> 
  mutate(percentage = round(syn_count/sum(syn_count), 3)*1) |>  
  filter(grepl("TuBu", pre_type)) |> 
  rename(er_id = post_pt_root_id,
         tubu_id = pre_pt_root_id,
         tubu_type = pre_type,
         er_input_percentage = percentage)|> 
  dplyr::select(er_type, er_id, tubu_id, tubu_type, er_input_percentage)  

tubu_df <- tubu_r_input_neuron |> 
  group_by(post_pt_root_id) |> 
  mutate(percentage = round(syn_count/sum(syn_count), 3)*1) |> 
  filter(grepl("MeTu", pre_type)) |>
  rename(tubu_id = post_pt_root_id,
         tubu_input_percentage = percentage) |> 
  dplyr::select( tubu_id,pre_pt_root_id, pre_type, tubu_input_percentage)  


###

tubu_tutu_df = tubu_r_input_neuron |> 
  group_by(post_pt_root_id) |> 
  mutate(percentage = round(syn_count/sum(syn_count), 3)*1) |> 
  filter(grepl("TuTu", pre_type)) |> 
  rename(tubu_id = post_pt_root_id,
         tutu_id = pre_pt_root_id,
         tutu_type = pre_type,
         tubu_input_percentage = percentage) |>  
  dplyr::select( tubu_id,tutu_id, tutu_type, tubu_input_percentage)


tutu_df <- tutu_input_neuron |>  
  group_by(post_pt_root_id) |>  
  mutate(percentage = round(syn_count/sum(syn_count), 3)*1) |>  
  filter(grepl("MeTu", pre_type)) |>
  rename(tutu_id = post_pt_root_id,
         tutu_input_percentage = percentage) |>  
  dplyr::select(tutu_id,pre_pt_root_id, pre_type, tutu_input_percentage) 

####

metur <- metu_r_input_me |> dplyr::select(post_pt_root_id, column, post_pt_position_x, post_pt_position_y, post_pt_position_z)
metur$hemisphere <- "r"

metur <- metur |> 
  mutate(all_id = paste(column,hemisphere, sep = "_"))

metul <-metu_l_input_me |>  dplyr::select(post_pt_root_id, column, post_pt_position_x, post_pt_position_y, post_pt_position_z)
metul$hemisphere <- "l"

metul <- metul |>   
  mutate(all_id = paste(column,hemisphere, sep = "_"))

metu_input_me <- rbind(metul,metur)

qw <- metu_input_me |>  dplyr::select(post_pt_root_id, all_id , column, hemisphere, column, post_pt_position_x, post_pt_position_y, post_pt_position_z)


vf_df <- er_df |>  left_join(tubu_df, by = "tubu_id") |>  
  left_join(qw, by = c("pre_pt_root_id"="post_pt_root_id")) |>  
  mutate(score = er_input_percentage * tubu_input_percentage) |>  
  group_by(er_id, column, hemisphere ) |>  
  mutate(value = sum(score)) |>  
  ungroup() |>  
  distinct(across(c(er_id, column, hemisphere)), .keep_all = T) |>  
  dplyr::select(er_type,er_id,column, all_id,hemisphere,value)



vf_df_minus <- er_df |>  left_join(tubu_tutu_df, by = "tubu_id") |> 
  left_join(tutu_df, by = c("tutu_id")) |> 
  left_join(qw, by = c("pre_pt_root_id"="post_pt_root_id")) |> 
  mutate(score = - er_input_percentage * tubu_input_percentage * tutu_input_percentage) |> 
  group_by(er_id, column, hemisphere ) |> 
  mutate(value = sum(score)) |> 
  ungroup() |> 
  distinct(across(c(er_id, column, hemisphere)), .keep_all = T) |> 
  dplyr::select(er_type,er_id,column, all_id,hemisphere,value)


# ER visual fields collection and population contours ---------------------


itt <- vf_df |> 
  distinct(er_id, .keep_all =T) |> 
  filter(er_type != "ER6") |> 
  pull(er_id) 

df_m = data.frame()
df_p = data.frame()

for (i in itt) {
  vf_df_sgn <- vf_df |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id, all_id,  value)
  
  
  vf_df_minus_sgn <- vf_df_minus |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id,all_id, value)
  
  
  er_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_sgn, by = "all_id") |> 
    mutate(n = round(coalesce(value, 0)*100)) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_minus_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_minus_sgn, by = "all_id") |> 
    mutate(n = abs(round(coalesce(value, 0)*100))) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_sgn_kde <- kde2d(x = er_sgn_hex_vf$phi,
                      y = er_sgn_hex_vf$theta,
                      h = c(60, 30),
                      n = c(180, 90), 
                      lims = c(-180, 180, -90, 90))
  
  er_sgn_kde_density_data <- expand.grid(x = er_sgn_kde$x, y = er_sgn_kde$y)
  er_sgn_kde_density_data$z <- as.vector(er_sgn_kde$z)
  er_sgn_kde_density_data$root_id = paste(i)
  er_sgn_kde_density_data$type = vf_df_sgn$er_type[1]
  
  if (nrow(er_minus_sgn_hex_vf) == 0){
    er_minus_sgn_kde_density_data <- expand.grid(x = seq(-180,180,by=2), y = seq(-90,90,by=2))
    er_minus_sgn_kde_density_data$z <- 0
    er_minus_sgn_kde_density_data$root_id = paste(i)
    er_minus_sgn_kde_density_data$type = vf_df_sgn$er_type[1]
  } else {
    er_minus_sgn_kde <- kde2d(x = er_minus_sgn_hex_vf$phi,
                              y = er_minus_sgn_hex_vf$theta,
                              h = c(60, 30),
                              n = c(180, 90), 
                              lims = c(-180, 180, -90, 90))
    er_minus_sgn_kde_density_data <- expand.grid(x = er_minus_sgn_kde$x, y = er_minus_sgn_kde$y)
    er_minus_sgn_kde_density_data$z <- as.vector(er_minus_sgn_kde$z)
    er_minus_sgn_kde_density_data$root_id = paste(i)
    er_minus_sgn_kde_density_data$type = vf_df_sgn$er_type[1]
  }
  df_p = rbind(df_p, er_sgn_kde_density_data)
  df_m = rbind(df_m, er_minus_sgn_kde_density_data)
}

max_type_dens_dir <- df_p |> 
  group_by(type) |> 
  summarise(max = max(z))

max_type_dens_indir <- df_m |> 
  group_by(type) |> 
  summarise(max = max(z)) |> 
  mutate(max = case_when(max == 0 ~ 1,
                         max != 0 ~ max))



for (i in itt) {
  
  er_type <- vf_df |> 
    filter(er_id == i) |>
    distinct(er_type) |>
    pull(er_type)
    
  
  vf_df_sgn <- vf_df |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id, all_id,  value)
  

  vf_df_minus_sgn <- vf_df_minus |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id,all_id, value)
  
  
  er_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_sgn, by = "all_id") |> 
    mutate(n = round(coalesce(value, 0)*100)) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_minus_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_minus_sgn, by = "all_id") |> 
    mutate(n = abs(round(coalesce(value, 0)*100))) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_sgn_kde <- kde2d(x = er_sgn_hex_vf$phi,
                      y = er_sgn_hex_vf$theta,
                      h = c(60, 30),
                      n = c(180, 90), 
                      lims = c(-180, 180, -90, 90))
  
  er_sgn_kde_density_data <- expand.grid(x = er_sgn_kde$x, y = er_sgn_kde$y)
  er_sgn_kde_density_data$z <- as.vector(er_sgn_kde$z)
  
  if (nrow(er_minus_sgn_hex_vf) == 0){
    er_minus_sgn_kde_density_data <- expand.grid(x = seq(-180,180,by=2), y = seq(-90,90,by=2))
    er_minus_sgn_kde_density_data$z <- 0
  } else {
    er_minus_sgn_kde <- kde2d(x = er_minus_sgn_hex_vf$phi,
                              y = er_minus_sgn_hex_vf$theta,
                              h = c(60, 30),
                              n = c(180, 90), 
                              lims = c(-180, 180, -90, 90))
    er_minus_sgn_kde_density_data <- expand.grid(x = er_minus_sgn_kde$x, y = er_minus_sgn_kde$y)
    er_minus_sgn_kde_density_data$z <- as.vector(er_minus_sgn_kde$z)
  }
  
  max_dir <- max_type_dens_dir|>
    filter(type == er_type) |> 
    pull(max)
  
  max_dir <- max_dir + max_dir*0.05
  
  plot_dir <- ggplot(er_sgn_kde_density_data, aes(x = x, y = y, z = z)) +
    geom_raster(aes(fill = z), interpolate = TRUE)+
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(
      colors = c("#FF000000", "#FF0000"),  
      guide = "colourbar",
      limits = c(0,max_dir)
      ) +
    labs(fill=paste(er_type, "density")) + 
    scale_y_continuous( limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    cowplot::theme_cowplot(font_size = 2,
                           line_size = 0.3,
                           font_family = "") +
    theme_void() 
  
  plot_dir_no_legend <- plot_dir +   
    theme(legend.position = "none")
  
  ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_red.png", sep = ""), plot = plot_dir_no_legend, width = 6, height = 3)
  
  legend_dir <- cowplot::get_legend(plot_dir)
  legend_dir <- ggdraw(legend_dir)
  
 # ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_direct_legend.pdf", sep = ""), plot = legend_dir, units = "mm", width = 30, height = 45)
  
  
  max_indir <- max_type_dens_indir|>
    filter(type == er_type) |> 
    pull(max)
  
  max_indir <- max_indir + max_indir*0.05
  
  plot_indir <- ggplot(er_minus_sgn_kde_density_data, aes(x = x, y = y, z = z)) +
    geom_raster(aes(fill = z), interpolate = TRUE)+
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(
      colors = c("#0000FF00", "#0000FF"), 
      guide = "colourbar",
      limits = c(0,max_indir)
      ) +
    labs(fill=paste(er_type, "density")) + 
    scale_y_continuous( limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    theme_void() 
  plot_indir_no_legend <- plot_indir + theme(legend.position = "none") 
  
  ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_blue.png", sep = ""), plot = plot_indir_no_legend, width = 6, height = 3)
  
  legend_indir <- cowplot::get_legend(plot_indir)
  legend_indir <- ggdraw(legend_indir)
  
  #ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_indirect_legend.pdf", sep = ""), plot = legend_indir, units = "mm", width = 30, height = 45)
  
  
  # Load the first image
  img1 <- readPNG(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_red.png", sep = ""))
  grob1 <- rasterGrob(img1, interpolate = TRUE)
  
  # Load the second image
  img2 <- readPNG(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_blue.png", sep = ""))
  grob2 <- rasterGrob(img2, interpolate = TRUE)
  
  
  ggplot() +
    annotation_custom(grob2, xmin = -180, xmax = 180, ymin = -90, ymax = 90) + 
    annotation_custom(grob1, xmin = -180, xmax = 180, ymin = -90, ymax = 90) +
    geom_point(data = eye_map_l, aes(phi, theta),  size = 0.17, stroke = 0, shape = 16, col = "grey25") +
    geom_point(data = eye_map_r, aes(phi, theta),  size = 0.17, stroke = 0, shape = 16,col = "grey25") +
    geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = 90), 
              fill = NA, color = "grey25", linewidth = 0.3)+
    xlab("azimuth") +
    ylab("elevation") +
    ggtitle(paste(er_type, itt, sep = "_"))+
    scale_y_continuous(breaks = seq(-90, 90, by = 45), limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    coord_fixed(1)+
    theme_minimal_grid(font_size = 6,
                           line_size = 0.17,
                           font_family = "") +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank())
  
  
  
  
  ggsave(paste("output/final/",vf_df_sgn$er_type[1],"_",i,"_vf.png", sep = "") , units = "mm", dpi = 1200, width =36.5, height = 18.25 )
  
}





itti <- vf_df |> 
  distinct(er_type, .keep_all =T) |> 
  filter(er_type != "ER6") |> 
  pull(er_type) 


for (i in itti) {
  azz = df_p |> filter(type == i)
  
  ggplot(azz, aes(x = x, y = y, col=root_id)) +
    geom_contour(aes(x = x, y = y, z = z, col=root_id),bins = 2,  size = 0.17) +
    geom_point(data = eye_map_l, aes(x = phi, y = theta), size = 0.17, stroke = 0, shape = 16, col = "grey25") +
    geom_point(data = eye_map_r, aes(x = phi, y = theta), size = 0.17, stroke = 0, shape = 16, col = "grey25") +
    geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = 90), 
              fill = NA, color = "grey25", linewidth = 0.3)+
    xlab("azimuth") +
    ylab("elevation") +
    ggtitle(paste(i))+
    
    scale_y_continuous(breaks = seq(-90, 90, by = 45), limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    coord_fixed(1)+
    theme_minimal_grid(font_size = 6,
                       line_size = 0.17,
                       font_family = "") +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          legend.position="none")
  
  ggsave(paste("output/final/",i,"_contour.png", sep = "") , units = "mm", dpi = 1200, width =36.5, height = 18.25)
  
}



# ER visual fields gallery ------------------------------------------------

for (i in itt) {
  
  er_type <- vf_df |> 
    filter(er_id == i) |>
    distinct(er_type) |>
    pull(er_type)
  
  
  vf_df_sgn <- vf_df |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id, all_id,  value)
  
  
  vf_df_minus_sgn <- vf_df_minus |> 
    filter(er_id == i) |> 
    dplyr::select(er_type,er_id,all_id, value)
  
  
  er_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_sgn, by = "all_id") |> 
    mutate(n = round(coalesce(value, 0)*100)) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_minus_sgn_hex_vf <- eye_map |> 
    left_join(vf_df_minus_sgn, by = "all_id") |> 
    mutate(n = abs(round(coalesce(value, 0)*100))) |> 
    filter(n > 0) |>        
    uncount(weights = n)|> 
    add_count(phi, theta) |> 
    mutate(dup = n > 1) |> 
    mutate(
      phi = if_else(dup, jitter(phi, amount = 2.5), phi),
      theta = if_else(dup, jitter(theta,amount = 2.5), theta)
    ) |> 
    dplyr::select(-n, -dup) |> 
    drop_na(theta)
  
  er_sgn_kde <- kde2d(x = er_sgn_hex_vf$phi,
                      y = er_sgn_hex_vf$theta,
                      h = c(60, 30),
                      n = c(180, 90), 
                      lims = c(-180, 180, -90, 90))
  
  er_sgn_kde_density_data <- expand.grid(x = er_sgn_kde$x, y = er_sgn_kde$y)
  er_sgn_kde_density_data$z <- as.vector(er_sgn_kde$z)
  
  if (nrow(er_minus_sgn_hex_vf) == 0){
    er_minus_sgn_kde_density_data <- expand.grid(x = seq(-180,180,by=2), y = seq(-90,90,by=2))
    er_minus_sgn_kde_density_data$z <- 0
  } else {
    er_minus_sgn_kde <- kde2d(x = er_minus_sgn_hex_vf$phi,
                              y = er_minus_sgn_hex_vf$theta,
                              h = c(60, 30),
                              n = c(180, 90), 
                              lims = c(-180, 180, -90, 90))
    er_minus_sgn_kde_density_data <- expand.grid(x = er_minus_sgn_kde$x, y = er_minus_sgn_kde$y)
    er_minus_sgn_kde_density_data$z <- as.vector(er_minus_sgn_kde$z)
  }
  
  max_dir <- max_type_dens_dir|>
    filter(type == er_type) |> 
    pull(max)
  
  max_dir <- max_dir + max_dir*0.05
  
  ggplot(er_sgn_kde_density_data, aes(x = x, y = y, z = z)) +
    geom_raster(aes(fill = z), interpolate = TRUE)+
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(
      colors = c("#FF000000", "#FF0000"),  
      guide = "colourbar",
      limits = c(0,max_dir)
    ) +
    scale_y_continuous( limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    cowplot::theme_cowplot(font_size = 2,
                           line_size = 0.3,
                           font_family = "") +
    theme_void() +
    theme(legend.position = "none")
  
  
  
  ggsave(paste("output/pngs/gallery/", vf_df_sgn$er_type[1],"_", i,"_red.png", sep = ""), width = 6, height = 3)
  
  #legend_dir <- cowplot::get_legend(plot_dir)
  #legend_dir <- ggdraw(legend_dir)
  
  # ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_direct_legend.pdf", sep = ""), plot = legend_dir, units = "mm", width = 30, height = 45)
  
  
  max_indir <- max_type_dens_indir|>
    filter(type == er_type) |> 
    pull(max)
  
  max_indir <- max_indir + max_indir*0.05
  
  ggplot(er_minus_sgn_kde_density_data, aes(x = x, y = y, z = z)) +
    geom_raster(aes(fill = z), interpolate = TRUE)+
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(
      colors = c("#0000FF00", "#0000FF"), 
      guide = "colourbar",
      limits = c(0,max_indir)
    ) +
    scale_y_continuous( limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    theme_void()+
    theme(legend.position = "none") 
  
  ggsave(paste("output/pngs/gallery/", vf_df_sgn$er_type[1],"_", i,"_blue.png", sep = ""), width = 6, height = 3)
  
  #legend_indir <- cowplot::get_legend(plot_indir)
  #legend_indir <- ggdraw(legend_indir)
  
  #ggsave(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_indirect_legend.pdf", sep = ""), plot = legend_indir, units = "mm", width = 30, height = 45)
  
  
  # Load the first image
  img1 <- readPNG(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_red.png", sep = ""))
  grob1 <- rasterGrob(img1, interpolate = TRUE)
  
  # Load the second image
  img2 <- readPNG(paste("output/pngs/", vf_df_sgn$er_type[1],"_", i,"_blue.png", sep = ""))
  grob2 <- rasterGrob(img2, interpolate = TRUE)
  
  
  ggplot() +
    geom_point(data = eye_map_l, aes(phi, theta), alpha = 0.8, size = 0.1, stroke = 0, shape = 16, col = "grey45") +
    geom_point(data = eye_map_r, aes(phi, theta), alpha = 0.8, size = 0.1, stroke = 0, shape = 16,col = "grey45") +
    annotation_custom(grob2, xmin = -180, xmax = 180, ymin = -90, ymax = 90) + 
    annotation_custom(grob1, xmin = -180, xmax = 180, ymin = -90, ymax = 90) +
    #geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = 90), 
    #          fill = NA, color = "grey25", linewidth = 0.3)+
    xlab("azimuth") +
    ylab("elevation") +
    ggtitle(paste(er_type, i, sep = "_"))+
    scale_y_continuous(breaks = seq(-90, 90, by = 45), limits = c(-90,90), expand = c(0,0))+
    scale_x_continuous(breaks = seq(180, -180, by = -45), limits = c(-180,180),expand = c(0,0))+
    coord_fixed(1)+
    theme_minimal_grid(line_size = 0.17,
                       font_family = "") +
    theme(text=element_text(size=4),
          axis.text=element_text(size=4),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.title = element_text(size = 4,
                                    hjust = 0.5, vjust = -2),
          axis.line = element_blank())
  # theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
  #       plot.title = element_text(size = 4,
  #                                 hjust = 0.5, vjust = -2),
  #       axis.title = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.line = element_blank(),
  #       axis.text = element_blank())
  
  
  
  
  ggsave(paste("output/final/gallery/",vf_df_sgn$er_type[1],"_",i,"_vf1.pdf", sep = "") , units = "mm", device = cairo_pdf, width =1.7*36.5, height = 1.7*18.25 )
  
}



pdf_folder <- "output/final/gallery/" # Change this to your folder path


# Define the names to group by
names_to_group <- c("ER2_ad", "ER2_b", "ER2_c", "ER3a_ad", "ER3d_a", "ER3d_b", "ER3d_c", "ER3d_d", "ER3m", "ER3p_ab", "ER3w_ab", "ER4d", "ER4m", "ER5", "ExR1")

# Function to extract the name from the file name
extract_name <- function(file_name) {
  name_parts <- str_split(file_name, "_")[[1]]
  
  # Combine all parts except the last two
  grouping_name <- str_c(name_parts[1:(length(name_parts) - 2)], collapse = "_")
  
  #name_parts <- unlist(strsplit(file_name, "_"))
  #grouping_name <- paste(name_parts[(length(name_parts)-2)], name_parts[(length(name_parts)-1)], sep = "_")
  if (grouping_name %in% names_to_group) {
    return(grouping_name)
  } else {
    return(NA)
  }
}




# Read all PDF files from the folder
pdf_files <- list.files(path = pdf_folder, pattern = "*.pdf", full.names = TRUE)
pdf_files_info <- data.frame(file = pdf_files, name = sapply(pdf_files, function(f) extract_name(basename(f))), stringsAsFactors = FALSE)

# Filter out files that do not match any name in names_to_group
pdf_files_info <- pdf_files_info[!is.na(pdf_files_info$name), ]

# Set PDF parameters
output_pdf_file <- "combined_output.pdf"
width <- 8.27 # A4 width in inches
height <- 11.69 # A4 height in inches
columns <- 3
rows <- 6 # Adjust based on how many images you want per page
dpi <- 600
# Calculate the size of each image cell
image_width <- width / columns
image_height <- height / rows

# Function to preserve aspect ratio
scale_image <- function(image, max_width, max_height) {
  img_width <- as.integer(image_info(image)$width)
  img_height <- as.integer(image_info(image)$height)
  
  scale <- min(max_width / img_width, max_height / img_height)
  
  image_resize(image, sprintf("%dx", as.integer(img_width * scale)))
}

# Function to start a new PDF page
start_new_page <- function() {
  grid.newpage()
}

# Create the PDF
pdf(output_pdf_file, width = width, height = height)

# Process each group
for (name in unique(pdf_files_info$name)) {
  cat("Processing name:", name, "\n")
  
  # Get files for the current name
  group_files <- pdf_files_info$file[pdf_files_info$name == name]
  
  # Read and process PDF pages for each file in the group
  group_images <- lapply(group_files, function(pdf) {
    pdf_pages <- pdf_info(pdf)$pages
    lapply(1:pdf_pages, function(page) {
      image_read(pdf_render_page(pdf, page = page, dpi = dpi))
    })
  })
  group_images <- unlist(group_images, recursive = FALSE)
  
  # Arrange images in the PDF
  for (i in seq_along(group_images)) {
    if ((i - 1) %% (columns * rows) == 0) {
      # Start a new page for the current name
      if (i != 1) {
        dev.off()
        pdf(output_pdf_file, width = width, height = height, onefile = TRUE)
      }
      start_new_page()
    }
    
    # Preserve aspect ratio while scaling image
    image <- scale_image(group_images[[i]], image_width * dpi, image_height * dpi) # High resolution
    
    # Calculate the position to center the image in its cell
    img_info <- image_info(image)
    x_pos <- ((i - 1) %% columns) * image_width + (image_width - (img_info$width / dpi)) / 2
    y_pos <- height - (((i - 1) %/% columns) %% rows + 1) * image_height + (image_height - (img_info$height / dpi)) / 2
    
    # Plot the image
    grid.raster(as.raster(image), x = unit(x_pos, "in"), y = unit(y_pos, "in"), width = unit(img_info$width / dpi, "in"), height = unit(img_info$height / dpi, "in"), just = c("left", "bottom"))
  }
}

# Close the last page
dev.off()





















