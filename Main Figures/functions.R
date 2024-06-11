## Collection of custom functions for working with flywire data.

# vec_norm() function -----------------------------------------------------

## Caluctates for a given vector (x) the vector magnitude and returns the 
## the normalized vector / unit vector. 

vec_norm <- function(x){
  vec_mag = apply(x, 1, function (x) sqrt(sum((x)^2)))
  x / vec_mag
}


# flywire_neuron_trans() function -----------------------------------------

## Performances a rigid transformation on a neuron mesh or a list of neuron meshes.
## Requires a center of rotation (center) and rotation axes (axis). 

flywire_neuron_trans <- function(neuron, axis, center){
  neuron_trans <- neuron
  for (i in 1:length(neuron)) {
    xyz <- neuron[[i]]$vb[1:3,]
    xyz_trans <- sweep(t(xyz), 2, center) %*% axis %>%  t()
    neuron_trans[[i]]$vb[1:3,] <- xyz_trans
  }
  return(neuron_trans)
}


# %nin%() function --------------------------------------------------------

## Helper function negating %in%.

`%nin%` = Negate(`%in%`)


# col_hex_grid() function -------------------------------------------------

## Transforming hexgrid matrix into dataframe.

col_hex_grid <- function(google_sheet){
  matrix = data.matrix(google_sheet)
  
  df = data.frame(col = as.numeric(),
                  x = as.numeric(),
                  y = as.numeric(),
                  v_pos = as.numeric(),
                  v_neg = as.numeric(),
                  p_pos = as.numeric(),
                  p_neg = as.numeric(),
                  q_pos = as.numeric(),
                  q_neg = as.numeric(),
                  h_pos = as.numeric(),
                  h_neg = as.numeric())
  
  cc = 0
  
  for (i in 3:(nrow(matrix)-2)) {
    for (j in 3:(ncol(matrix)-2)) {
      df[cc+j-1, "col"] = matrix[i,j]
      df[cc+j-1, "x"] = j
      df[cc+j-1, "y"] = -i
      df[cc+j-1, "v_pos"] = matrix[i-2,j]
      df[cc+j-1, "v_neg"] = matrix[i+2,j]
      df[cc+j-1, "p_pos"] = matrix[i-1,j-1]
      df[cc+j-1, "p_neg"] = matrix[i+1,j+1]
      df[cc+j-1, "q_pos"] = matrix[i-1,j+1]
      df[cc+j-1, "q_neg"] = matrix[i+1,j-1]
      df[cc+j-1, "h_pos"] = matrix[i,j+2]
      df[cc+j-1, "h_neg"] = matrix[i,j-2]
    }
    cc =  cc + length(2:(ncol(matrix)-1)) 
  }
  
  df = df %>% filter(!is.na(col)) 
  return(df)
}


# approx3d() function -----------------------------------------------------


approx3d <- function(x1, y1, z1, x2, y2, z2, v) {
  vvv = c()
  for (i in 1:length(v)) {
    vv = c((c(x1, y1, z1) - c(x2, y2, z2)) * v[i] + c(x2, y2, z2)) 
    ss = matrix(vv, ncol = 3, byrow = F)
    vvv = rbind(vvv, ss)
  }
  return(vvv)
}


# neuron_class_syn_table() function ---------------------------------------

## processes the "project_synapse_table" table into output and input synapse tables 
## with named input and output neurons based on the garner_kind_labels. Returns three
## dataframes in a list.


neuron_syn_table <- function(class, hemisphere = c("L", "R"), filter_column = "class") {
  
  # Ensure filter_column is either "class" or "label"
  if (!filter_column %in% c("class", "label")) {
    stop("filter_column must be either 'class' or 'label'")
  }
  
  filter_column_sym <- sym(filter_column)
  
  class_id <- garner_kind_labels |> 
    filter(!!filter_column_sym == {{class}} & hemisphere %in% {{hemisphere}})
  
  dynamic_col_name <- tolower(paste0(class, "_type"))
  
  class_input <- project_synapse_table |> 
    filter(post_pt_root_id %in% class_id$root_id) |> 
    left_join(class_id[, c("label", "root_id")], by = c("post_pt_root_id" = "root_id")) |> 
    rename(!!dynamic_col_name := label) |> 
    left_join(project_labels, by = c("pre_pt_root_id" = "root_id")) |> 
    rename(pre_type = label) 
  
  class_output <- project_synapse_table |> 
    filter(pre_pt_root_id %in% class_id$root_id) |> 
    left_join(class_id[, c("label", "root_id")], by = c("pre_pt_root_id" = "root_id")) |> 
    rename(!!dynamic_col_name := label) |> 
    left_join(project_labels, by = c("pre_pt_root_id" = "root_id")) |> 
    rename(pre_type = label) 
  
  return_list <- list(
    class_id = class_id,
    class_input = class_input,
    class_output = class_output
  )
  
  if (length(hemisphere) == 1) {
    names(return_list) <- tolower(c(
      paste(class, hemisphere, "id", sep = "_"),
      paste(class, hemisphere, "input", sep = "_"),
      paste(class, hemisphere, "output", sep = "_")
    ))
  } else {
    names(return_list) <- tolower(c(
      paste(class, "id", sep = "_"),
      paste(class, "input", sep = "_"),
      paste(class, "output", sep = "_")
    ))
  }
  
  return(return_list)
}


# conn_matrix_plot() function ---------------------------------------------

## generates a pdf file in the output folder with connectivity matrices (heatmaps)

conn_matrix_plot <- function(metutype = c("MeTu1", "MeTu2a", "MeTu2b",
                                          "MeTu3a", "MeTu3b", "MeTu3c",
                                          "MeTu4a", "MeTu4b", "MeTu4c",
                                          "MeTu4d"),
                             posttype = c("TuBu01", "TuBu02", "TuBu03",
                                          "TuBu04", "TuBu05", "TuBu06",
                                          "TuBu07", "TuBu08", "TuBu09",
                                          "TuBu10"),
                             cluster_n =1,
                             population = F,
                             file_name = "MeTu_dend_matrix"){
  
  top_input <- unlist(ls_proofed[metutype]) |> 
    setdiff("Unidentified")
  
  if(isTRUE(population)){
    df_me <- metu_r_me_input_type
    df_aotu <- metu_r_aotu_output_type
  } else{
    df_me <- metu_r_me_input_type_proofed
    df_aotu <- metu_r_aotu_output_type_proofed
  }
  
  metu_r_input_cluster <- df_me |>
    filter(metu_type %in%  metutype & pre_type %in% top_input) |> 
    dplyr::select(post_pt_root_id, metu_type, pre_type, percentage) |> 
    pivot_wider(names_from = pre_type, values_from = percentage, values_fill = 0) |> 
    rename(root_id = post_pt_root_id) |> 
    ungroup()
  
  metu_r_output_cluster <- df_aotu |> 
    filter(metu_type %in% metutype & post_type %in% posttype ) |>  
    dplyr::select(pre_pt_root_id, metu_type, post_type, percentage) |> 
    pivot_wider(names_from = post_type, values_from = percentage, values_fill = 0) |> 
    rename(root_id = pre_pt_root_id) |> 
    ungroup()
  
  #metu_cluster <- cbind(metu_r_input_cluster, metu_r_output_cluster[,-c(1:2)] )
  metu_cluster <- metu_r_input_cluster |> 
    dplyr::select(-metu_type) |> 
    inner_join(metu_r_output_cluster, by = "root_id") 
  
  metu_r_input_cluster <- metu_r_input_cluster |> 
    filter(root_id %in% metu_cluster$root_id) 
  
  metu_r_output_cluster <- metu_r_output_cluster |> 
    filter(root_id %in% metu_cluster$root_id)
  
  metu_cluster_mat <- metu_cluster |> 
    dplyr::select(-c(root_id,metu_type)) |>  
    as.data.frame()
  rownames(metu_cluster_mat) <- metu_cluster$root_id
  
  metu_row_dend <- as.dendrogram(hclust(dist(metu_cluster_mat)))
  metu_row_dend <- color_branches(metu_row_dend, k = cluster_n) 
  
  metu_row_labels <- rownames(metu_cluster_mat)
  metu_row_labels <- metu_r |>
    filter(root_id %in% metu_row_labels) |> 
    pull(label, root_id)
  
  
  metu_in_cluster_mat <- metu_r_input_cluster |> 
    dplyr::select(-c(root_id,metu_type))  |> 
    as.data.frame()
  rownames(metu_in_cluster_mat) <- metu_r_input_cluster$root_id
  
  metu_in_row_dend <- as.dendrogram(hclust(dist(metu_in_cluster_mat)))
  metu_in_row_dend <- color_branches(metu_in_row_dend, k = cluster_n) 
  
  metu_in_row_labels <- rownames(metu_in_cluster_mat)
  metu_in_row_labels  <- metu_r |>
    pull(label, root_id)
  
  
  metu_out_cluster_mat <- metu_r_output_cluster |>
    dplyr::select(-c(root_id,metu_type))  |>
    as.data.frame()
  rownames(metu_out_cluster_mat) <- metu_r_output_cluster$root_id
  
  metu_out_row_dend <- as.dendrogram(hclust(dist(metu_out_cluster_mat)))
  metu_out_row_dend <- color_branches(metu_out_row_dend, k = cluster_n) 
  
  metu_out_row_labels <- rownames(metu_out_cluster_mat)
  metu_out_row_labels <- metu_r |>
    filter(root_id %in% metu_out_row_labels) |> 
    pull(label, root_id)
  
  #row_split <- cutree(hclust(dist(metu_cluster_mat)), k = cluster_n)
  if(cluster_n >1){
    row_split <- cluster_n
  }else{
    row_split <- NULL
  }
  
  if(length(metutype) > 6 & isTRUE(population)){
    pdf_height <- 40
  } else if(length(metutype) < 6 & isTRUE(population)){
    pdf_height <- 20
  } else if(length(metutype) > 6 & !isTRUE(population)){
    pdf_height <- 10
  } else{
    pdf_height <- 5
  }
  
  if(length(metutype) > 6){
    metu_row_labels <- sapply(metu_row_labels, function(x) {
      sub("MeTu", "", x)
    })
  } else {
    metu_row_labels <- sapply(metu_row_labels, function(x) {
      gsub("MeTu\\d+", "", x)
    })
  }
  
  mycols_in <- colorRamp2(breaks = c(0, max(metu_in_cluster_mat)), colors = c("white", "darkcyan"))
  mycols_out <- colorRamp2(breaks = c(0, max(metu_out_cluster_mat)), colors = c("white", "darkred"))
  
  # Define the dimensions for square tiles
  tile_size <- unit(2, "mm")
  height_in <- nrow(metu_in_cluster_mat) * tile_size
  width_in <- ncol(metu_in_cluster_mat) * tile_size
  height_out <- nrow(metu_out_cluster_mat) * tile_size
  width_out <- ncol(metu_out_cluster_mat) * tile_size
  
  # Create the heatmaps with smaller font size for column labels and square tiles
  hm_metu_in <- Heatmap(metu_in_cluster_mat,
                        show_column_dend = FALSE,
                        cluster_rows = metu_row_dend,
                        col = mycols_in,
                        row_dend_reorder = FALSE,
                        row_dend_width = unit(1, "cm"),
                        width = width_in,
                        height = height_in,
                        row_labels = metu_row_labels[rownames(metu_cluster_mat)],
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 6),
                        show_heatmap_legend = TRUE,
                        row_split = row_split,  
                        gap = unit(2, "mm"))
  
  hm_metu_out <- Heatmap(metu_out_cluster_mat,
                         show_column_dend = FALSE,
                         cluster_rows = metu_row_dend,
                         col = mycols_out,
                         row_dend_reorder = FALSE,
                         row_dend_width = unit(1, "cm"),
                         width = width_out,
                         height = height_out,
                         row_labels = metu_row_labels[rownames(metu_cluster_mat)],
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         show_heatmap_legend = TRUE,
                         row_split = row_split,
                         gap = unit(2, "mm"))
  
  # Save the combined heatmap to a PDF file with adjusted layout
  pdf(file = paste("output/", file_name,".pdf"), width = 6, height = pdf_height)
  
  # Combine and draw the heatmaps, adjusting the layout
  draw(hm_metu_out + hm_metu_in, heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE, gap = unit(2, "mm"))
  
  dev.off()
  
  # return(row_split)
}
