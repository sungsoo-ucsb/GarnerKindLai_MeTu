# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 19:49:49 2022

@author: Dustin Garner
"""


import os
import time
import math
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from fafbseg import flywire
import flycode.utils as utils
import flycode.readfiles as readfiles
import flycode.mapping as mapping
import flycode.neuprint_reading as neuread
import flycode.reduction as reduction
import flycode.specific as specific
import flycode.maintenance as maintenance
import flycode.meshes as meshes
import flycode.proofreading as proofreading
import flycode.synapse_density as density
import flycode.metu_comparison as comparison
import flycode.annotations as annotations
import flycode.volumes as volumes
import flycode.flywire_functions as fw
import flycode.figures as figures


line_width = 0.2#0.0645
plt.rcParams["font.family"] = ["Arial", "sans-serif"]
plt.rcParams["axes.linewidth"] = line_width
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["xtick.major.width"] = line_width
plt.rcParams["xtick.minor.width"] = line_width
plt.rcParams["ytick.major.width"] = line_width
plt.rcParams["ytick.minor.width"] = line_width


#Easy reference type lists.
all_mt_types = ["1","2a","2b","3a","3b","3c","4a","4b","4c","4d"]
mt1_l = ["MeTu1_L"]
mt1_r = ["MeTu1_R"]
mt2_l = [f"MeTu2{x}_L" for x in "ab"]
mt2_r = [f"MeTu2{x}_R" for x in "ab"]
mt3_l = [f"MeTu3{x}_L" for x in "abc"]
mt3_r = [f"MeTu3{x}_R" for x in "abc"]
mt4_l = [f"MeTu4{x}_L" for x in "abcd"]
mt4_r = [f"MeTu4{x}_R" for x in "abcd"]
mt_l = mt1_l + mt2_l + mt3_l + mt4_l
mt_r = mt1_r + mt2_r + mt3_r + mt4_r
all_mt_l = mt_l
all_mt_r = mt_r + ["MeTu_incomplete_R"]
all_mt = mt_l + mt_r
mt_l_general = [f"MeTu{x}_L" for x in "1234"]
mt_r_general = [f"MeTu{x}_R" for x in "1234"]
tb_pl_l = ["TuBu08_L"]
tb_pl_r = ["TuBu08_R"]
pc_tb_types = "16"
tb_pc_l = [f"TuBu0{x}_L" for x in pc_tb_types]
tb_pc_r = [f"TuBu0{x}_R" for x in pc_tb_types]
a_tb_types = ["07","09","10"]
tb_a_l = [f"TuBu{x}_L" for x in a_tb_types]
tb_a_r = [f"TuBu{x}_R" for x in a_tb_types]
m_tb_types = "2345"
tb_m_l = [f"TuBu0{x}_L" for x in m_tb_types]
tb_m_r = [f"TuBu0{x}_R" for x in m_tb_types]
all_tb_types = [f"0{x!s}" if x<10 else str(x) for x in range(1,11)]
tb_l = [f"TuBu{x}_L" for x in all_tb_types]
tb_r = [f"TuBu{x}_R" for x in all_tb_types]
er_pl_l = ["ER4d_L"]
er_pl_r = ["ER4d_R"]
pc_ring_types = ["4m","5"]
er_pc_l = [f"ER{x}_L" for x in pc_ring_types]
er_pc_r = [f"ER{x}_R" for x in pc_ring_types]
a_ring_types = ["3w_ab","2_ad","2_b","2_c"]
er_a_l = [f"ER{x}_L" for x in a_ring_types] 
er_a_r = [f"ER{x}_R" for x in a_ring_types]
er_a = [f"ER{x}_{y}" for x in a_ring_types for y in "LR"]
m_ring_types = ["3a_ad","3m","3d_a","3d_b","3d_c","3d_d","3p_ab"]
er_m_l = [f"ER{x}_L" for x in m_ring_types]
er_m_r = [f"ER{x}_R" for x in m_ring_types]
all_ring_types = ["2_ad","2_b","2_c",
                  "3a_ad","3d_a","3d_b","3d_c","3d_d","3m","3p_ab","3w_ab",
                  "4d","4m",
                  "5"]
er_l = [f"ER{x}_L" for x in all_ring_types]
er_r = [f"ER{x}_R" for x in all_ring_types]
er_all = [f"ER{x}_{y}" for x in all_ring_types for y in "LR"]
retino_ring_order = ["4d"]+pc_ring_types+a_ring_types+m_ring_types
er_l_retino = [f"ER{x}_L" for x in retino_ring_order]
er_r_retino = [f"ER{x}_R" for x in retino_ring_order]

aotu46 = aotu046 = [f"AOTU046_{x}" for x in "LR"]
aotu46_indvid = [f"AOTU046_{x}_{y}" for x in "LR" for y in "12"]
tutu = [f"TuTuB_{x}_{y}" for x in "ab" for y in "LR"]
avp_all = all_mt_l+all_mt_r+tb_l+["TuBu_misc_L"]+tb_r+er_l+er_r+aotu46+tutu


#Easy Access Comparison Types
gen_metu = [f"MeTu{x}" for x in ["1","2a","2b","3ab","3c","4a","4b","4c","4d",
                                 "4e","4f"]]
gen_tubu = [f"TuBu{x}" for x in all_tb_types]
gen_ring = [f"ER{x}" for x in all_ring_types]
gen_bihem = ["AOTU046", "TuTuB_a", "TuTuB_b"]


upstream_neurons = ['Dm2_R','DmDRA1_R','Mi15_R','R7_R','R7DRA_R','TmY31_R',
                    'cM01c_R','cM02b_R','cM08a_R','cM12_L','cM13_L','MTe01a_R'] +\
                   [f'{x}_{y}' for y in "LR" for x in ['MeMeDRA','aMe19b',
                                                       'MeMe_e01','MeMe_e02',
                                                       'MeMe_e10',]] 
    

upstream_with_mi1 = upstream_neurons[:13] + ["Mi1_L", "Mi1_R"] + \
            upstream_neurons[13:]


epg = [f"EPG_R{x!s}" for x in range(8,0,-1)]+\
      [f"EPG_L{x!s}" for x in range(1,9)]
epg_retino = [f"EPG_{x}{y}" for x,y in zip("LR"*8, 
        np.array([[a,b] for a,b in zip(range(1,9),range(8,0,-1))]).flatten())]
er_all_plus_lal = ["ER1_L", "ER1_R"] + er_all[:8] + ["ER3a_bc_L", "ER3a_bc_R"] +\
            er_all[8:] + ["ER6_L", "ER6_R"]
er_all_all = er_all_plus_lal[:12] + ["ER3_misc_R"] + er_all_plus_lal[12:]

exr = [f"ExR{x!s}_{y}" for x in range(1,9) for y in "LR"]

epg = [f"EPG_R{x!s}" for x in range(8,0,-1)]+\
      [f"EPG_L{x!s}" for x in range(1,9)]
epg_retino = [f"EPG_{x}{y}" for x,y in zip("LR"*8, 
        np.array([[a,b] for a,b in zip(range(1,9),range(8,0,-1))]).flatten())]

pen1 = [f"PEN1_R{x!s}" for x in range(9,1,-1)]+\
       [f"PEN1_L{x!s}" for x in range(2,10)]
pen2 = [f"PEN2_R{x!s}" for x in range(9,1,-1)]+\
       [f"PEN2_L{x!s}" for x in range(2,10)]
peg = [f"PEG_R{x!s}" for x in range(9,0,-1)]+\
      [f"PEG_L{x!s}" for x in range(1,10)]
el = [f"EL_{x}" for x in "LR"]
epgt = [f"EPGt_{x}" for x in "LR"]

eb = epg + pen1 + pen2 + peg + el + epgt
eb_broad = [f"{x}_{y}" for x in ["EPG","PEN1","PEN2","PEG","EL","EPGt"] 
                       for y in "LR"]

all_in_paper = [f"MeTu{x}_{y}" for x in all_mt_types for y in "LR"] +\
               ["MeTu_incomplete_R"] +\
               [f"TuBu{x}_{y}" for x in all_tb_types for y in "LR"] +\
               ["TuBu_misc_L"] +\
               er_all_all +\
               tutu +\
               aotu46 +\
               upstream_with_mi1 +\
               exr +\
               eb
        
mt1_r_outlier_coords = np.array([[202927, 71601, 4917]])
mt2_r_outlier_coords = np.array([[170645, 42981, 1977]])
mt4_r_outlier_coords = np.array([[202103, 56303, 5479], 
                                    [203472, 61664, 5513],
                                    [201806, 62921, 4839],])

#ring_outlier_coords = np.array([[118158, 44346, 2494]])
fig_size_mt_tb = (1.5, 1.5)
bihem_size = (3, 3)
fig_size_mt1 = (1.3, 1.3)
fig_size_mt_tb_pl = (3.6, 3.6)
fig_size_tb_er_pl_r = (1.3, 1.3)
fig_size_mt2 = (1.1, 1.1)
fig_size_mt_tb_pc = (3.5, 3.5)
fig_size_tb_er_pc = (1.3, 1.3)
fig_size_mt3 = (1.1, 1.1)
fig_size_mt_tb_a = (3.5, 3.5)
fig_size_tb_er_a = (1.3, 1.3)
fig_size_mt4 = (1.1, 1.1)
fig_size_mt_tb_m = (3.45, 3.45)
fig_size_tb_er_m = (2.2, 2.2)

regions = {"AME_L": 1, "LO_L": 2, "NO": 3, "BU_L": 4, "PB": 5, "LH_L": 6,
           "LAL_L": 7, "SAD": 8, "CAN_L": 9, "AMMC_L": 10, "ICL_L": 11,
           "VES_L": 12, "IB_L": 13, "ATL_L": 14, "CRE_L": 15, "MB_PED_L": 16,
           "MB_VL_L": 17, "MB_ML_L": 18, "FLA_L": 19, "LOP_L": 20, "EB": 21,
           "AL_L": 22, "ME_L": 23, "FB": 24, "SLP_L": 25, "SIP_L": 26,
           "SMP_L": 27, "AVLP_L": 28, "PVLP_L": 29, "WED_L": 30, "PLP_L": 31,
           "AOTU_L": 32, "GOR_L": 33, "MB_CA_L": 34, "SPS_L": 35, "IPS_L": 36,
           "SCL_L": 37, "EPA_L": 38, "GNG": 39, "PRW": 40, "GA_L": 41,
           "AME_R": 42, "LO_R": 43, "BU_R": 44, "LH_R": 45, "LAL_R": 46,
           "CAN_R": 47, "AMMC_R": 48, "ICL_R": 49, "VES_R": 50, "IB_R": 51,
           "ATL_R": 52, "CRE_R": 53, "MB_PED_R": 54, "MB_VL_R": 55, 
           "MB_ML_R": 56, "FLA_R": 57, "LOP_R": 58, "AL_R": 59, "ME_R": 60,
           "SLP_R": 61, "SIP_R": 62, "SMP_R": 63, "AVLP_R": 64, "PVLP_R": 65,
           "WED_R": 66, "PLP_R": 67, "AOTU_R": 68, "GOR_R": 69, "MB_CA_R": 70,
           "SPS_R": 71, "IPS_R": 72, "SCL_R": 73, "EPA_R": 74, "GA_R": 75}
for i in regions.copy():
    regions[regions[i]] = i

plot_folder = ""
directory = os.path.dirname(__file__)
directory = os.path.join(directory, "flycode")
for i in ["Excel-Plots", "Generated-Figures", "Importable-Coords",
          "Meshes", "Readable"]:
    subfolders = [x.path for x in os.scandir(directory) if x.is_dir()]
    temp_directory = os.path.join(directory, i)
    if not temp_directory in subfolders:
        os.mkdir(temp_directory)




"""
# -------- #
# Figure 1 #
# -------- #
plot_folder = "Fig 1"

#Fig. 1c
colors = ["TPurple", "TPurple", "TBlue", "TRed", "TGreen", "TOrange2", 
          "TOrange2", "TRed", "TGreen", "TPink"]
for i, j in zip(maps, colors):
    density.plot_maps(maps[i], color=j, plot_name=f"{i} Blurrier", blur=10, 
                      max_value=0.157)

#Fig. 1di
mt_tb_r = mapping.ConnectionMap(mt_r_general, tb_r, "AOTU_R")
mt_tb_r.make_type_plots(plot_name = "MeTu_R to TuBu_R", 
                        plot_folder=plot_folder, fig_size=fig_size_mt_tb)

#Fig. 1dii
mt_tb_l = mapping.ConnectionMap(mt_l_general, tb_l, "AOTU_L")
mt_tb_l.make_type_plots(plot_name = "MeTu_L to TuBu_L",
                        plot_folder=plot_folder, fig_size=fig_size_mt_tb)

#Fig. 1diii
mt_tb_r = mapping.ConnectionMap(mt_r, tb_r, "AOTU_R") #Specific
mt_tb_r.make_type_plots(plot_name="Subtyped MeTu_R to TuBu_R",
                        plot_folder=plot_folder)

#Fig. 1div
mt_tb_l = mapping.ConnectionMap(mt_l, tb_l, "AOTU_L") #Specific
mt_tb_l.make_type_plots(plot_name="Subtyped MeTu_L to TuBu_L",
                        plot_folder=plot_folder)


# --------- #
# Figure S1 #
# --------- #
plot_folder = "Fig S1"

#Fig. S1biii and S1biv
f1_scores = volumes.compare_f1_scores()
volumes.plot_f1(f1_scores, plot_folder=plot_folder)
volumes.get_all_t_tests(f1_scores) #Gets t-test p-scores.


# --------- #
# Figure S2 #
# --------- #
#Figure this one out


# --------- #
# Figure S3 #
# --------- #
plot_folder = "Fig S3"

#Fig. S3a
bihem_aotu_r = mapping.ConnectionMap(aotu46+tutu+mt_r, aotu46+tutu+mt_r+tb_r, "AOTU_R")
bihem_aotu_r.make_type_plots(\
            plot_name="Bihemispheric Connections (AOTU_R)", plot_folder=plot_folder,
            fig_size=bihem_size)

#Fig. S3b
bihem_aotu_l = mapping.ConnectionMap(aotu46+tutu+mt_l, aotu46+tutu+mt_l+tb_l, "AOTU_L")
bihem_aotu_l.make_type_plots(\
            plot_name="Bihemispheric Connections (AOTU_L)", plot_folder=plot_folder,
            fig_size=bihem_size)

#Fig. S3c
bihem_bu_r = mapping.ConnectionMap(aotu46+tb_r, aotu46+tb_r+er_r, "BU_R")
bihem_bu_r.make_type_plots(\
            plot_name="Bihemispheric Connections (BU_R)", plot_folder=plot_folder,
            fig_size=bihem_size)

#Fig. S3d
bihem_bu_l = mapping.ConnectionMap(aotu46+tb_l, aotu46+tb_l+er_l, "BU_L")
bihem_bu_l.make_type_plots(\
            plot_name="Bihemispheric Connections (BU_L)", plot_folder=plot_folder,
            fig_size=bihem_size)


# --------- #
# Figure S4 #
# --------- #
plot_folder = "Fig S4"

#Fig. S4
eb_inter = mapping.ConnectionMap(er_all_plus_lal+exr+eb_broad, 
                                 er_all_plus_lal+exr+eb_broad, "EB")
eb_inter.make_type_plots(plot_name="EB Interconnectivity",
                         plot_folder=plot_folder, fig_size=(8, 8))


# --------- #
# Figure S5 #
# --------- #
plot_folder = "Fig S5"

#Fig. S5ai
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_me = mapping.ConnectionMap(mt1_r, mt1_r, "ME_R", exclude=mt1_r_outliers)
mt1_r_inter_me.make_connectivity_plots(plot_name="MT1_R Interconnectivity (ME_R)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt1)

#Fig. S5aii
mt1_l_inter_me = mapping.ConnectionMap(mt1_l, mt1_l, "ME_L")
mt1_l_inter_me.make_connectivity_plots(plot_name="MT1_L Interconnectivity (ME_L)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt1)

#Fig. S5bi
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_aotu = mapping.ConnectionMap(mt1_r, mt1_r, "AOTU_R",)
mt1_r_inter_aotu.make_connectivity_plots(plot_name="MT1_R Interconnectivity (AOTU_R)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt1)

#Fig. S5bii
mt1_l_inter_aotu = mapping.ConnectionMap(mt1_l, mt1_l, "AOTU_L")
mt1_l_inter_aotu.make_connectivity_plots(plot_name="MT1_L Interconnectivity (AOTU_L)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt1)

#Fig. S5ci
mt_tb_pl_r = mapping.ConnectionMap(mt1_r, tb_pl_r, "AOTU_R")
mt_tb_pl_r.make_connectivity_plots(plot_name="MeTu1_R to TuBu_R",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S5cii
mt_tb_pl_l = mapping.ConnectionMap(mt1_l, tb_pl_l, "AOTU_L")
mt_tb_pl_l.make_connectivity_plots(plot_name="MeTu1_L to TuBu_L",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S5di
tb_er_pl_r = mapping.ConnectionMap(tb_pl_r, er_pl_r, "BU_R")
tb_er_pl_r.make_connectivity_plots(plot_name="TuBu_PL_R to ER_PL_R",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_tb_er_pl_r)

#Fig. S5dii
fig_size_tb_er_pl_l = (1.15, 1.15)
tb_er_pl_l = mapping.ConnectionMap(tb_pl_l, er_pl_l, "BU_L")
tb_er_pl_l.make_connectivity_plots(plot_name="TuBu_PL_L to ER_PL_L",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_tb_er_pl_l)


# --------- #
# Figure S6 #
# --------- #
plot_folder = "Fig S6"

#Fig. S6a
line_widths = specific.get_bihem_weight_line_width()
for i in line_widths:
    print(i, line_widths[i])

#Fig. S6b, S6di-ii
specific.make_bihem_pie_charts(plot_folder=plot_folder)

#Fig. S6c
comparison.nt_by_types(["AOTU046_L_1", "AOTU046_L_2", "AOTU046_R_1", "AOTU046_R_2"],
        palette=["#e6194B", "#EBC400", "#3E7748", "#4400dd"],
        plot_names=["AOTU046 NT Predictions"], plot_folder=plot_folder,
        fig_size=(1, 0.75))

#Fig. S6ei
comparison.nt_by_types(["TuTuB_a_L", "TuTuB_a_R"], 
                       plot_names=["TuTuB_a NT Predictions"], 
                       plot_folder=plot_folder,
                       palette=["#e6194B", "#3E7748"])

#Fig. S6eii
comparison.nt_by_types(["TuTuB_b_L", "TuTuB_b_R"],
                       plot_names=["TuTuB_b NT Predictions"], 
                       plot_folder=plot_folder,
                       palette=["#e6194B", "#3E7748"])

#Fig. S6f
specific.tutu_comparison(plot_name="TuTu Synapse Counts Comparison", 
                         plot_folder=plot_folder)


# --------- #
# Figure S7 #
# --------- #
plot_folder = "Fig S7"

#Fig. S7ai
mt2_r_outliers = fw.locs_to_segments(mt2_r_outlier_coords)
mt2_r_inter_me = mapping.ConnectionMap(mt2_r, mt2_r, "ME_R", exclude = mt2_r_outliers)
mt2_r_inter_me.make_connectivity_plots(plot_name="MT2_R Interconnectivity (ME_R)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt2)

#Fig. S7aii
mt2_l_inter_me = mapping.ConnectionMap(mt2_l, mt2_l, "ME_L")
mt2_l_inter_me.make_connectivity_plots(plot_name="MT2_L Interconnectivity (ME_L)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt2)

#Fig. S7bi
mt2_r_inter_aotu = mapping.ConnectionMap(mt2_r, mt2_r, "AOTU_R",)
mt2_r_inter_aotu.make_connectivity_plots(plot_name="MT2_R Interconnectivity (AOTU_R)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt2)

#Fig. S7bii
mt2_l_inter_aotu = mapping.ConnectionMap(mt2_l, mt2_l, "AOTU_L")
mt2_l_inter_aotu.make_connectivity_plots(plot_name="MT2_L Interconnectivity (AOTU_L)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt2)

#Fig. S7ci
mt_tb_pc_r = mapping.ConnectionMap(mt2_r, tb_pc_r, "AOTU_R")
mt_tb_pc_r.make_connectivity_plots(plot_name="MeTu2_R to TuBu_R",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S7cii
mt_tb_pc_l = mapping.ConnectionMap(mt2_l, tb_pc_l, "AOTU_L")
mt_tb_pc_l.make_connectivity_plots(plot_name="MeTu2_L to TuBu_L",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S7di
tb_er_pc_r = mapping.ConnectionMap(tb_pc_r, er_pc_r, "BU_R")
tb_er_pc_r.make_connectivity_plots(plot_name="TuBu_PC_R to ER_PC_R",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_tb_er_pc)

#Fig. S7dii
tb_er_pc_l = mapping.ConnectionMap(tb_pc_l, er_pc_l, "BU_L")
tb_er_pc_l.make_connectivity_plots(plot_name="TuBu_PC_L to ER_PC_L",
                                   plot_folder=plot_folder,
                                   fig_size=fig_size_tb_er_pc)


# --------- #
# Figure S8 #
# --------- #
plot_folder = "Fig S8"

#Fig. S8ai
mt3_r_inter_me = mapping.ConnectionMap(mt3_r, mt3_r, "ME_R")
mt3_r_inter_me.make_connectivity_plots(plot_name="MT3_R Interconnectivity (ME_R)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt3)

#Fig. S8aii
mt3_l_inter_me = mapping.ConnectionMap(mt3_l, mt3_l, "ME_L")
mt3_l_inter_me.make_connectivity_plots(plot_name="MT3_L Interconnectivity (ME_L)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt3)

#Fig. S8bi
mt3_r_inter_aotu = mapping.ConnectionMap(mt3_r, mt3_r, "AOTU_R")
mt3_r_inter_aotu.make_connectivity_plots(plot_name="MT3_R Interconnectivity (AOTU_R)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt3)

#Fig. S8bii
mt3_l_inter_aotu = mapping.ConnectionMap(mt3_l, mt3_l, "AOTU_L")
mt3_l_inter_aotu.make_connectivity_plots(plot_name="MT3_L Interconnectivity (AOTU_L)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt3)

#Fig. S8ci
mt_tb_a_r = mapping.ConnectionMap(mt3_r, tb_a_r, "AOTU_R")
mt_tb_a_r.make_connectivity_plots(plot_name="MeTu3_R to TuBu_R",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_mt_tb_a)

#Fig. S8cii
mt_tb_a_l = mapping.ConnectionMap(mt3_l, tb_a_l, "AOTU_L")
mt_tb_a_l.make_connectivity_plots(plot_name="MeTu3_L to TuBu_L",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_mt_tb_a)

#Fig. S8di
tb_er_a_r = mapping.ConnectionMap(tb_a_r, er_a_r, "BU_R")
tb_er_a_r.make_connectivity_plots(plot_name="uBu_A_R to ER_A_R",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_tb_er_a)

#Fig. S8dii
tb_er_a_l = mapping.ConnectionMap(tb_a_l, er_a_l, "BU_L")
tb_er_a_l.make_connectivity_plots(plot_name="TuBu_A_L to ER_A_L",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_tb_er_a)

#Fig. S8e
specific.sm17_map_by_dv_axis(plot_name= 
        "MeTu3bc Vertical Positioning by Connections to Sm17",
        plot_folder=plot_folder)

#Fig. S8g
specific.mt3_pre_connections(["Sm23", "MeMeDRA", "Mi15"], 
                             plot_name="MeTu3 Presynaptic", plot_folder=plot_folder)


# --------- #
# Figure S9 #
# --------- #
plot_folder = "Fig S9"

#Fig. S9ai
mt4_r_outliers = fw.locs_to_segments(mt4_r_outlier_coords)
mt4_r_inter_me = mapping.ConnectionMap(mt4_r, mt4_r, "ME_R", exclude = mt4_r_outliers)
mt4_r_inter_me.make_connectivity_plots(plot_name="MT4_R Interconnectivity (ME_R)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt4)

#Fig. S9aii
mt4_l_inter_me = mapping.ConnectionMap(mt4_l, mt4_l, "ME_L",)
mt4_l_inter_me.make_connectivity_plots(plot_name="MT4_L Interconnectivity (ME_L)",
                                       plot_folder=plot_folder,
                                       fig_size=fig_size_mt4)

#Fig. S9bi
mt4_r_inter_aotu = mapping.ConnectionMap(mt4_r, mt4_r, "AOTU_R",)
mt4_r_inter_aotu.make_connectivity_plots(plot_name="MT4_R Interconnectivity (AOTU_R)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt4)

#Fig. S9bii
mt4_l_inter_aotu = mapping.ConnectionMap(mt4_l, mt4_l, "AOTU_L",)
mt4_l_inter_aotu.make_connectivity_plots(plot_name="MT4_L Interconnectivity (AOTU_L)",
                                         plot_folder=plot_folder,
                                         fig_size=fig_size_mt4)

#Fig. S9ci
mt_tb_m_r = mapping.ConnectionMap(mt4_r, tb_m_r, "AOTU_R")
mt_tb_m_r.make_connectivity_plots(plot_name="MeTu4_R to TuBu_R",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_mt_tb_m)

#Fig. S9cii
mt_tb_m_l = mapping.ConnectionMap(mt4_l, tb_m_l, "AOTU_L")
mt_tb_m_l.make_connectivity_plots(plot_name="MeTu4_L to TuBu_L",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_mt_tb_m)

#Fig. S9di
tb_er_m_r = mapping.ConnectionMap(tb_m_r, er_m_r, "BU_R")
tb_er_m_r.make_connectivity_plots(plot_name="TuBu_M_R to ER_M_R",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_tb_er_m)

#Fig. S9dii
tb_m_l.insert(1, "TuBu_misc_L")
tb_er_m_l = mapping.ConnectionMap(tb_m_l, er_m_l, "BU_L")
tb_er_m_l.make_connectivity_plots(plot_name="TuBu_M_L to ER_M_L",
                                  plot_folder=plot_folder,
                                  fig_size=fig_size_tb_er_m)
tb_m_l.remove("TuBu_misc_L")


# ---------- #
# Figure S12 #
# ---------- #
plot_folder = "Fig S12"

#Fig. S12a
neuread.plot_comparison(gen_metu,
        datasets={"FAFB": "LR", "Hemibrain": "R", "FIB-SEM": "R"},
        plot_name="MeTu Neuron Counts Comparison", 
        plot_folder=plot_folder,
        y_ticks=np.arange(0, 140, 20))

#Fig. S12b
neuread.plot_comparison(gen_tubu, 
        plot_name="TuBu Neuron Counts Comparison", 
        plot_folder=plot_folder,
        y_ticks=np.arange(0, 14, 2))

#Fig. S12c
neuread.plot_comparison(gen_ring,
        datasets={"FAFB": "LR", "Hemibrain": "LR"},
        plot_name="Ring Neuron Counts Comparison", 
        plot_folder=plot_folder,
        y_ticks=np.arange(0, 14, 2), fig_size = (2.4, 1.25))

#Fig. S12d
neuread.plot_comparison(gen_bihem, 
        datasets={"FAFB": "LR", "Hemibrain": "LR"},
        plot_name="Bihemispheric Neuron Counts Comparison", 
        plot_folder=plot_folder,
        y_ticks=np.arange(0, 3, 1))

#Fig. S12e
neuread.ratio_plot(neuread.Ratio.METU_TO_TUBU, plot_name=\
    "Ratio of MeTu Neurons to TuBu Neurons per Hemisphere",
    plot_folder=plot_folder,
    y_ticks=np.arange(0, 18, 3))

#Fig. S12f
neuread.ratio_plot(neuread.Ratio.TUBU_TO_RING, plot_name=\
    "Ratio of TuBu Neurons to Ring Neurons per Hemisphere",
    plot_folder=plot_folder,
    y_ticks=np.arange(0,5,1))

#Fig. S12g
neuread.lobula_counts(just_metu4=True,
        plot_name="Lobula Synapse Count in Hemibrain MeTu4",
        plot_folder=plot_folder,)

#Fig. S12h
specific.full_comparison([f"MeTu4{x}" for x in "abcd"], "LO", 
        plot_name="Lobula Synapse Count per MeTu4 Neuron", 
        plot_folder=plot_folder,
        y_axis="Lobula Synapse Count")

#Fig. S12i_i-ii
comparison.mt4e_partner_comparison("Partner Type Comparison", 
                                   plot_folder=plot_folder,
                                   save_figure=True)

#Fig. S12i_iii-viii
comparison.mt4e_dorsal_comparison("Comparing",
                                  plot_folder=plot_folder,
                                  save_figure=True)


# ---------- #
# Figure S13 #
# ---------- #
plot_folder = "Fig S13"

#Fig. S13ai
comparison.compare_metu("Medulla Presynapse Count", plot_folder=plot_folder)

#Fig. S13aii
comparison.compare_metu("Medulla Postsynapse Count", plot_folder=plot_folder)

#Fig. S13aiii
comparison.compare_metu("Column Count", plot_folder=plot_folder)

#Fig. S13aiv
comparison.compare_metu("AOTU Presynapse Count", plot_folder=plot_folder)

#Fig. S13av
comparison.compare_metu("AOTU Postsynapse Count", plot_folder=plot_folder)

#Fig. S13avi
comparison.compare_metu("Ellipse Ratio", plot_folder=plot_folder)

#Fig. S13bi-viii
for i in ["Medulla Postsynapse Count", "AOTU Postsynapse Count",
          "Medulla Presynapse Count", "AOTU Presynapse Count",
          "Column Count", "Ellipse Ratio", "Ellipse Major Axis Length (nm)",
          "Ellipse Minor Axis Length (nm)"]:
    for j in ["A-P", "V-D"]:
        comparison.scatter_plots(x=f"Relative Offset from Medulla Centroid ({j})",
                                 y=i, plot_name=f"{i} {j}", plot_folder=plot_folder)


# ---------- #
# Figure S14 #
# ---------- #
plot_folder = "Fig S14"

#Fig. S14a-j
metu_nts = comparison.nt_by_types(\
        [f"{x}_R" for x in comparison.colors], 
        palette=list(comparison.colors.values()), 
        plot_names=[f"{x} Neurotransmitter Prediction" for x in comparison.colors.keys()],
        plot_folder=plot_folder,
        save_figure=True,
        fig_size=(1,0.75), 
        separate_plots=True)


# ---------- #
# Figure S15 #
# ---------- #
plot_folder = "Fig S15"


#Fig. S15bi
mapping.add_broad_type("TuBu_R", tb_r)
mapping.add_broad_type("TuTuB_a", [f"TuTuB_a_{x}" for x in "LR"])
mapping.add_broad_type("TuTuB_b", [f"TuTuB_b_{x}" for x in "LR"])
mapping.add_broad_type("AOTU046", [f"AOTU046_{x}" for x in "LR"])
all_tubu_weights = mapping.get_total_weight(mt_r+["TuTuB_a", "TuTuB_b", "AOTU046"],
                                             ["TuBu_R"], region="AOTU_R")

#Fig. S15bii
mapping.add_broad_type("ER_R", er_r)
mapping.add_broad_type("AOTU046", [f"AOTU046_{x}" for x in "LR"])
all_ring_weights = mapping.get_total_weight(tb_r+["AOTU046"], ["ER_R"], region="BU_R")

#Fig. S15c
ring_to_epg = mapping.ConnectionMap(er_all, epg, "EB")
ring_to_epg.make_type_plots("Ring to EPG", plot_folder=plot_folder, fig_size=(3.0, 4.0))


# ------------ #
# Spreadsheets #
# ------------ #

#Download "Labels" from https://codex.flywire.ai/api/download and store it
#in "flycode/Readable" as "Codex Labels 783.csv"
label_df = proofreading.all_codex_names(all_in_paper)
utils.write_excel(label_df, "Codex Naming History")

percent_df = proofreading.full_percent_spread(all_in_paper)
utils.write_excel(percent_df, "FlyWire Consortium Edit Record")


# ------ #
# Colors #
# ------ #
#"#FF0000" - Red for presynaptic
#"#00FFFF" - Cyan for postsynaptic

#"#D62728" - Red for compare figure and MeTu1
#"#1F77B4" - Blue for compare figure and MeTu2
#"#2CA02C" - Green for compare figure and MeTu3
#"#FF8F0E" - Better Yellow for compare figure and MeTu4
#"#FF7F0E" - Orange
#"#BC27D6" - Purple
#"#D627B5" - Pink
#"#27B0D6" - Cyan


"""







