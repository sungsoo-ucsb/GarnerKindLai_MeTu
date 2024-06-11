<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 19:49:49 2022

@author: dusti
"""


import os
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
import flycode.graphics as graphics


desktop = "/Users/dustin/Desktop/"

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

aotu46 = [f"AOTU046_{x}" for x in "LR"]
aotu46_indvid = [f"AOTU046_{x}_{y}" for x in "LR" for y in "12"]
tutu = [f"TuTuB_{x}_{y}" for x in "ab" for y in "LR"]
avp_all = all_mt_l+all_mt_r+tb_l+["TuBu_misc_L"]+tb_r+er_l+er_r+aotu46+tutu


#Easy Access Comparison Types
gen_metu = [f"MeTu{x}" for x in ["1","2a","2b","3ab","3c","4a","4b","4c","4d"]]
gen_tubu = [f"TuBu{x}" for x in all_tb_types]
gen_ring = [f"ER{x}" for x in all_ring_types]
gen_bihem = ["AOTU046", "TuTuB_a", "TuTuB_b"]

upstream_neurons = ['Dm2_R','DmDRA1_R','MeLoPlp_unknown_1_R','MeLo_unknown_1_R',
    'MeMeDRA_R','MeMeDRA_L','MeMe_unknown_1_L','MeMe_unknown_1_R','MeMe_unknown_2_L',
    'MeMe_unknown_2_R','MeMe_unknown_3_L','MeMe_unknown_3_R','MeSps_unknown_1_R',
    'Mi15_R','Mti_DRA2_R',
    'Mti_unknown_1_R','Mti_unknown_2_R','Mti_unknown_3_R','Mti_unknown_4_R',
    'Mti_unknown_5_R','Mti_unknown_6_R',
    'Mti_unknown_7_R','Mti_unknown_8_R','Mti_unknown_9_R','Mti_unknown_10_R',
    "R7_R", "R7DRA_R",
    'SpsMe_unknown_1_R','SpsMe_unknown_2_R','SpsSpsMe_unknown_1_R']

upstream_with_mi1 = upstream_neurons[:13] + ["Mi1_L", "Mi1_R"] + \
            upstream_neurons[13:]


epg = [f"EPG_R{x!s}" for x in range(8,0,-1)]+\
      [f"EPG_L{x!s}" for x in range(1,9)]
epg_retino = [f"EPG_{x}{y}" for x,y in zip("LR"*8, 
        np.array([[a,b] for a,b in zip(range(1,9),range(8,0,-1))]).flatten())]
er_all_plus_lal = ["ER1_L", "ER1_R"] + er_all[:8] + ["ER3a_bc_L", "ER3a_bc_R"] +\
            er_all[8:]
er_all_all = er_all_plus_lal[:12] + ["ER3_misc_R"] + er_all_plus_lal[12:]

exr = [f"ExR{x!s}_{y}" for x in range(1,9) for y in "LR"]
exr = ["ExR1_L", "ExR1_R", "ExR2_L", "ExR2_R", "ExR3_L", "ExR3_R", "ExR4_L", "ExR4_R",
       "ExR5_L", "ExR5_R", "ExR6_L", "ExR6_R", "ExR7_L", "ExR7_R", "ExR8_L", "ExR8_R"]

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
                                 [201806, 62921, 4839], 
                                 [200665, 61972, 4648],
                                 [202760, 59834, 5178]])
mt4_l_outlier_coords = np.array([[58831, 63835, 6098], 
                                 [60314, 53680, 6226]])
ring_outlier_coords = np.array([[118158, 44346, 2494]])
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



"""
--------
Figure 1
--------
#Fig. 1c
#figure this one out

#Fig. 1di
mt_tb_r = mapping.ConnectionMap(mt_r_general, tb_r, "AOTU_R")
mt_tb_r.make_type_plots(plot_name = "Fig 1/MeTu_R to TuBu_R",
                        fig_size = fig_size_mt_tb)

#Fig. 1dii
mt_tb_l = mapping.ConnectionMap(mt_l_general, tb_l, "AOTU_L")
mt_tb_l.make_type_plots(plot_name = "Fig 1/MeTu_L to TuBu_L",
                        fig_size = fig_size_mt_tb)

#Fig. 1diii
mt_tb_r = mapping.ConnectionMap(mt_r, tb_r, "AOTU_R") #Specific
mt_tb_r.make_type_plots(plot_name = "Fig 1/Subtyped MeTu_R to TuBu_R")

#Fig. 1div
mt_tb_l = mapping.ConnectionMap(mt_l, tb_l, "AOTU_L") #Specific
mt_tb_l.make_type_plots(plot_name = "Fig 1/Subtyped MeTu_L to TuBu_L")


---------
Figure S1
---------
#Fig. S1b
volumes.plot_f1(volumes.compare_f1_scores()) #Get this to work

#Fig. S1ei
bihem_aotu_r = mapping.ConnectionMap(aotu46+tutu+mt_r, aotu46+tutu+mt_r+tb_r, "AOTU_R")
bihem_aotu_r.make_type_plots(\
            plot_name="Fig S1/Type Plots/Bihemispheric Connections (AOTU_R)",
            fig_size=bihem_size)

#Fig. S1eii
bihem_aotu_l = mapping.ConnectionMap(aotu46+tutu+mt_l, aotu46+tutu+mt_l+tb_l, "AOTU_L")
bihem_aotu_l.make_type_plots(\
            plot_name="Fig S1/Type Plots/Bihemispheric Connections (AOTU_L)",
            fig_size=bihem_size)

#Fig. S1fi
bihem_bu_r = mapping.ConnectionMap(aotu46+tb_r, aotu46+tb_r+er_r, "BU_R")
bihem_bu_r.make_type_plots(\
            plot_name="Fig S1/Type Plots/Bihemispheric Connections (BU_R)",
            fig_size=bihem_size)

#Fig. S1fii
bihem_bu_l = mapping.ConnectionMap(aotu46+tb_l, aotu46+tb_l+er_l, "BU_L")
bihem_bu_l.make_type_plots(\
            plot_name="Fig S1/Type Plots/Bihemispheric Connections (BU_L)",
            fig_size=bihem_size)

#Fig. S1g
ring_inter = mapping.ConnectionMap(er_all, er_all, "EB")
ring_inter.make_type_plots(
            plot_name="Fig S1/Type Plots/Ring Neuron Interconnectivity",
            fig_size=(3.2, 3.2))

#Fig. S1h
neuread.plot_comparison(gen_metu, 
        plot_name="Fig S1/Neuprint Comparison/MeTu Neuron Counts", 
        y_ticks=np.arange(0, 140, 20))

#Fig. S1i
neuread.plot_comparison(gen_tubu, 
        plot_name="Fig S1/Neuprint Comparison/TuBu Neuron Counts", 
        y_ticks=np.arange(0, 14, 2))

#Fig. S1j
neuread.plot_comparison(gen_ring, 
        plot_name="Fig S1/Neuprint Comparison/Ring Neuron Counts", 
        y_ticks=np.arange(0, 14, 2), fig_size = (2.4, 1.25))

#Fig. S1k
neuread.plot_comparison(gen_bihem, 
        plot_name="Fig S1/Neuprint Comparison/Bihemispheric Neuron Counts", 
        y_ticks=np.arange(0, 3, 1))

#Fig. S1l
neuread.ratio_plot(neuread.Ratio.METU_TO_TUBU, plot_name=\
    "Fig S1/Neuprint Comparison/Ratio of MeTu Neurons to TuBu Neurons per Hemisphere",\
    y_ticks=np.arange(0, 18, 3))

#Fig. S1m
neuread.ratio_plot(neuread.Ratio.TUBU_TO_RING, plot_name=\
    "Fig S1/Neuprint Comparison/Ratio of TuBu Neurons to Ring Neurons per Hemisphere",\
    y_ticks=np.arange(0,5,1))

    
---------
Figure S2
---------
#Fig. S2ai
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_me = mapping.ConnectionMap(mt1_r, mt1_r, "ME_R", exclude=mt1_r_outliers)
mt1_r_inter_me.make_connectivity_plots(plot_name="Fig S2/MT1_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt1)

#Fig. S2aii
mt1_l_inter_me = mapping.ConnectionMap(mt1_l, mt1_l, "ME_L")
mt1_l_inter_me.make_connectivity_plots(plot_name="Fig S2/MT1_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt1)

#Fig. S2bi
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_aotu = mapping.ConnectionMap(mt1_r, mt1_r, "AOTU_R", exclude=mt1_r_outliers)
mt1_r_inter_aotu.make_connectivity_plots(plot_name="Fig S2/MT1_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt1)

#Fig. S2bii
mt1_l_inter_aotu = mapping.ConnectionMap(mt1_l, mt1_l, "AOTU_L")
mt1_l_inter_aotu.make_connectivity_plots(plot_name="Fig S2/MT1_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt1)

#Fig. S2ci
mt_tb_pl_r = mapping.ConnectionMap(mt1_r, tb_pl_r, "AOTU_R")
mt_tb_pl_r.make_connectivity_plots(plot_name="Fig S2/MeTu1_R to TuBu_R",
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S2cii
mt_tb_pl_l = mapping.ConnectionMap(mt1_l, tb_pl_l, "AOTU_L")
mt_tb_pl_l.make_connectivity_plots(plot_name="Fig S2/MeTu1_L to TuBu_L",
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S2di
tb_er_pl_r = mapping.ConnectionMap(tb_pl_r, er_pl_r, "BU_R")
tb_er_pl_r.make_connectivity_plots(plot_name="Fig S2/TuBu_PL_R to ER_PL_R",
                                   fig_size=fig_size_tb_er_pl_r)

#Fig. S2dii
fig_size_tb_er_pl_l = (1.15, 1.15)
tb_er_pl_l = mapping.ConnectionMap(tb_pl_l, er_pl_l, "BU_L")
tb_er_pl_l.make_connectivity_plots(plot_name="Fig S2/TuBu_PL_L to ER_PL_L",
                                   fig_size=fig_size_tb_er_pl_l)

#Fig. S2h
specific.make_bihem_pie_charts()

#Fig. S2i
comparison.nt_by_types(["AOTU046_L_1", "AOTU046_L_2", "AOTU046_R_1", "AOTU046_R_2"],
        palette=["#e6194B", "#EBC400", "#3E7748", "#4400dd"],
        plot_names=["Fig S2/AOTU046 NT Predictions"], fig_size=(1, 0.75))



---------
Figure S3
---------
#Fig. S3ai
mt2_r_outliers = fw.locs_to_segments(mt2_r_outlier_coords)
mt2_r_inter_me = mapping.ConnectionMap(mt2_r, mt2_r, "ME_R", exclude = mt2_r_outliers)
mt2_r_inter_me.make_connectivity_plots(plot_name="Fig S3/MT2_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt2)

#Fig. S3aii
mt2_l_inter_me = mapping.ConnectionMap(mt2_l, mt2_l, "ME_L")
mt2_l_inter_me.make_connectivity_plots(plot_name="Fig S3/MT2_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt2)

#Fig. S3bi
mt2_r_inter_aotu = mapping.ConnectionMap(mt2_r, mt2_r, "AOTU_R", exclude = mt2_r_outliers)
mt2_r_inter_aotu.make_connectivity_plots(plot_name="Fig S3/MT2_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt2)

#Fig. S3bii
mt2_l_inter_aotu = mapping.ConnectionMap(mt2_l, mt2_l, "AOTU_L")
mt2_l_inter_aotu.make_connectivity_plots(plot_name="Fig S3/MT2_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt2)

#Fig. S3ci
mt_tb_pc_r = mapping.ConnectionMap(mt2_r, tb_pc_r, "AOTU_R")
mt_tb_pc_r.make_connectivity_plots(plot_name="Fig S3/MeTu2_R to TuBu_R",
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S3cii
mt_tb_pc_l = mapping.ConnectionMap(mt2_l, tb_pc_l, "AOTU_L")
mt_tb_pc_l.make_connectivity_plots(plot_name="Fig S3/MeTu2_L to TuBu_L",
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S3di
tb_er_pc_r = mapping.ConnectionMap(tb_pc_r, er_pc_r, "BU_R")
tb_er_pc_r.make_connectivity_plots(plot_name="Fig S3/TuBu_PC_R to ER_PC_R",
                                   fig_size=fig_size_tb_er_pc)

#Fig. S3dii
tb_er_pc_l = mapping.ConnectionMap(tb_pc_l, er_pc_l, "BU_L")
tb_er_pc_l.make_connectivity_plots(plot_name="Fig S3/TuBu_PC_L to ER_PC_L",
                                   fig_size=fig_size_tb_er_pc)

#Fig. S3hi
comparison.nt_by_types(["TuTuB_a_L", "TuTuB_a_R"], 
                       plot_names=["Fig S3/TuTuB_a NT Predictions"],
                       palette=["#e6194B", "#3E7748"])

#Fig. S3hii
comparison.nt_by_types(["TuTuB_b_L", "TuTuB_b_R"],
                       plot_names=["Fig S3/TuTuB_b NT Predictions"],
                       palette=["#e6194B", "#3E7748"])



---------
Figure S4
---------
#Fig. S4ai
mt3_r_inter_me = mapping.ConnectionMap(mt3_r, mt3_r, "ME_R")
mt3_r_inter_me.make_connectivity_plots(plot_name="Fig S4/MT3_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt3)

#Fig. S4aii
mt3_l_inter_me = mapping.ConnectionMap(mt3_l, mt3_l, "ME_L")
mt3_l_inter_me.make_connectivity_plots(plot_name="Fig S4/MT3_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt3)

#Fig. S4bi
mt3_r_inter_aotu = mapping.ConnectionMap(mt3_r, mt3_r, "AOTU_R")
mt3_r_inter_aotu.make_connectivity_plots(plot_name="Fig S4/MT3_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt3)

#Fig. S4bii
mt3_l_inter_aotu = mapping.ConnectionMap(mt3_l, mt3_l, "AOTU_L")
mt3_l_inter_aotu.make_connectivity_plots(plot_name="Fig S4/MT3_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt3)

#Fig. S4ci
mt_tb_a_r = mapping.ConnectionMap(mt3_r, tb_a_r, "AOTU_R")
mt_tb_a_r.make_connectivity_plots(plot_name="Fig S4/MeTu3_R to TuBu_R",
                                  fig_size=fig_size_mt_tb_a)

#Fig. S4cii
mt_tb_a_l = mapping.ConnectionMap(mt3_l, tb_a_l, "AOTU_L")
mt_tb_a_l.make_connectivity_plots(plot_name="Fig S4/MeTu3_L to TuBu_L",
                                  fig_size=fig_size_mt_tb_a)

#Fig. S4di
tb_er_a_r = mapping.ConnectionMap(tb_a_r, er_a_r, "BU_R")
tb_er_a_r.make_connectivity_plots(plot_name="Fig S4/TuBu_A_R to ER_A_R",
                                  fig_size=fig_size_tb_er_a)

#Fig. S4dii
tb_er_a_l = mapping.ConnectionMap(tb_a_l, er_a_l, "BU_L")
tb_er_a_l.make_connectivity_plots(plot_name="Fig S4/TuBu_A_L to ER_A_L",
                                  fig_size=fig_size_tb_er_a)

#Fig. S4e
specific.mti1_map_by_dv_axis(plot_name= 
        "Fig S4/MeTu3bc Vertical Positioning by Connections to Mti_unknown_1")

#Fig. S4g
specific.mt3_pre_connections(["Mti_DRA2", "MeMeDRA", "Mi15"], 
                             plot_name="Fig S4/MeTu3 Presynaptic")



---------
Figure S5
---------
#Fig. S5ai
mt4_r_outliers = fw.locs_to_segments(mt4_r_outlier_coords)
mt4_r_inter_me = mapping.ConnectionMap(mt4_r, mt4_r, "ME_R", exclude = mt4_r_outliers)
mt4_r_inter_me.make_connectivity_plots(plot_name="Fig S5/MT4_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt4)

#Fig. S5aii
mt4_l_outliers = fw.locs_to_segments(mt4_l_outlier_coords)
mt4_l_inter_me = mapping.ConnectionMap(mt4_l, mt4_l, "ME_L", exclude = mt4_l_outliers)
mt4_l_inter_me.make_connectivity_plots(plot_name="Fig S5/MT4_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt4)

#Fig. S5bi
mt4_r_inter_aotu = mapping.ConnectionMap(mt4_r, mt4_r, "AOTU_R", exclude = mt4_r_outliers)
mt4_r_inter_aotu.make_connectivity_plots(plot_name="Fig S5/MT4_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt4)

#Fig. S5bii
mt4_l_inter_aotu = mapping.ConnectionMap(mt4_l, mt4_l, "AOTU_L", exclude = mt4_l_outliers)
mt4_l_inter_aotu.make_connectivity_plots(plot_name="Fig S5/MT4_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt4)

#Fig. S5ci
mt_tb_m_r = mapping.ConnectionMap(mt4_r, tb_m_r, "AOTU_R")
mt_tb_m_r.make_connectivity_plots(plot_name="Fig S5/MeTu4_R to TuBu_R",
                                  fig_size=fig_size_mt_tb_m)

#Fig. S5cii
mt_tb_m_l = mapping.ConnectionMap(mt4_l, tb_m_l, "AOTU_L")
mt_tb_m_l.make_connectivity_plots(plot_name="Fig S5/MeTu4_L to TuBu_L",
                                  fig_size=fig_size_mt_tb_m)

#Fig. S5di
tb_er_m_r = mapping.ConnectionMap(tb_m_r, er_m_r, "BU_R")
tb_er_m_r.make_connectivity_plots(plot_name="Fig S5/TuBu_M_R to ER_M_R",
                                  fig_size=fig_size_tb_er_m)

#Fig. S5dii
tb_m_l.insert(1, "TuBu_misc_L")
tb_er_m_l = mapping.ConnectionMap(tb_m_l, er_m_l, "BU_L")
tb_er_m_l.make_connectivity_plots(plot_name="Fig S5/TuBu_M_L to ER_M_L",
                                  fig_size=fig_size_tb_er_m)
tb_m_l.remove("TuBu_misc_L")

#Fig. S5g
neuread.lobula_counts(just_metu4=True,\
        plot_name="Fig S5/Lobula Synapse Count in Hemibrain MeTu4")

#Fig. S5h
specific.full_comparison([f"MeTu4{x}" for x in "abcd"], "LO", 
        plot_name="Fig S5/Lobula Synapse Count per MeTu4 Neuron", 
        y_axis="Lobula Synapse Count")
    

---------
Figure S6
---------
#Fig. S6ai
comparison.compare_metu("Medulla Presynapse Counts")

#Fig. S6aii
comparison.compare_metu("Medulla Postsynapse Counts")

#Fig. S6aiii
comparison.compare_metu("Column Counts")

#Fig. S6aiv
comparison.compare_metu("AOTU Presynapse Counts")

#Fig. S6av
comparison.compare_metu("AOTU Postsynapse Counts")

#Fig. S6avi
comparison.compare_metu("Ellipse Ratio")

#Fig. S6bi
comparison.scatter_plots(x="Position Along Medulla Ventral-Dorsal Axis",\
   y="Medulla Presynapse Counts", plot_folder = "ME Pre Counts to DV/")
comparison.scatter_plots(x="Position Along Medulla Posterior-Anterior Axis",\
   y="Medulla Presynapse Counts", plot_folder = "ME Pre Counts to AP/")

#Fig. S6bii
comparison.scatter_plots(x="Position Along Medulla Ventral-Dorsal Axis",\
   y="Medulla Postsynapse Counts", plot_folder="ME Post Counts to DV/")
comparison.scatter_plots(x="Position Along Medulla Posterior-Anterior Axis",\
   y="Medulla Postsynapse Counts", plot_folder="ME Post Counts to AP/")

#Fig. S6biii
comparison.scatter_plots(x="Position Along Medulla Ventral-Dorsal Axis",\
   y="Column Counts", plot_folder="Column Counts to DV/")

#Fig. S6biv
comparison.scatter_plots(x="Position Along Medulla Ventral-Dorsal Axis",\
   y="AOTU Presynapse Counts", plot_folder="AOTU Pre Counts to DV/")
comparison.scatter_plots(x="Position Along Medulla Posterior-Anterior Axis",\
   y="AOTU Presynapse Counts", plot_folder="AOTU Pre Counts to AP/")

#Fig. S6bv
comparison.scatter_plots(x="Position Along Medulla Ventral-Dorsal Axis",\
   y="AOTU Postsynapse Counts", plot_folder="AOTU Post Counts to DV/")
comparison.scatter_plots(x = "Position Along Medulla Posterior-Anterior Axis",\
   y="AOTU Postsynapse Counts", plot_folder="AOTU Post Counts to AP/")

#Fig. S6bvi
comparison.scatter_plots(x="Position Along Medulla Posterior-Anterior Axis",\
   y="Column Counts", plot_folder="Column Counts to AP/")

#Fig. S6c
comparison.nt_by_types(
    [f"{x}_R" for x in comparison.colors], 
    palette=list(comparison.colors.values()), 
    plot_names=[f"Fig S6/Neurotransmitters/{x}" for x in comparison.colors],
    save_figure=True,
    fig_size=(1,0.75), 
    separate_plots=True)


---------
Figure S7
---------
#Fig. S7
ring_outliers = fw.locs_to_segments(ring_outlier_coords)
eb_inter = mapping.ConnectionMap(er_all_plus_lal+exr+eb_broad, 
                                 er_all_plus_lal+exr+eb_broad, "EB")
eb_inter.make_type_plots(plot_name = "Fig S7/EB Interconnectivity",
                         fig_size = (10, 10))



------------
Spreadsheets
------------

Download "Labels" from https://codex.flywire.ai/api/download before running:
proofreading.all_codex_names(all_in_paper)

"""




"""
============
Unused Plots
============

-------------
Synaptic Maps
-------------
#Adjust font_size to ? for these
tb_er_r_alpha = mapping.ConnectionMap(tb_r, er_r, "BU_R") #Alphabetical
tb_er_r_alpha.make_connectivity_plots(\
        plot_name="TuBu to Ring/Alphabetical TuBu_R to ER_R")

    
tb_l.insert(2, "TuBu_misc_L")
tb_er_l_alpha = mapping.ConnectionMap(tb_l, er_l, "BU_L") #Alphabetical
tb_er_l_alpha.make_connectivity_plots(\
        plot_name="TuBu to Ring/Alphabetical TuBu_L to ER_L")
tb_l.remove("TuBu_misc_L")


tb_er_r_retino = mapping.ConnectionMap(tb_r, er_r_retino, "BU_R") #Retinotopical
tb_er_r_retino.make_connectivity_plots(\
        plot_name="TuBu to Ring/Retinotopic TuBu_R to ER_R")

    
tb_l.insert(2, "TuBu_misc_L")
tb_er_l_retino = mapping.ConnectionMap(tb_l, er_l_retino, "BU_L") #Retinotopical
tb_er_l_retino.make_connectivity_plots(\
        plot_name="TuBu to Ring/Retinotopic TuBu_L to ER_L")
tb_l.remove("TuBu_misc_L")


er_a_inter = mapping.ConnectionMap(er_a, er_a, "EB")
er_a_inter.make_connectivity_plots(\
        plot_name="Ring Interconnectivity/Anterior Ring Interconnectivity")

    
ring_outliers = fw.locs_to_segments(ring_outlier_coords)
er_m_inter = mapping.ConnectionMap(er_m_l+er_m_r, er_m_l+er_m_r, "EB", 
                                   exclude = ring_outliers)
er_m_inter.make_connectivity_plots(\
        plot_name="Ring Interconnectivity/Medial Ring Interconnectivity")
#er_m_inter.to_excel(file_name = "Medial Ring Interconnectivity")





---------
Type Maps
---------
tb_er_r = mapping.ConnectionMap(tb_r, er_r, "BU_R") #Alphabetical
tb_er_r.make_type_plots(plot_name = "Type Connections/Alphabetical TuBu_R to ER_R")


tb_er_l = mapping.ConnectionMap(tb_l, er_l, "BU_L") #Alphabetical
tb_er_l.make_type_plots(plot_name = "Type Connections/Alphabetical TuBu_L to ER_L")


tb_er_r_retino = mapping.ConnectionMap(tb_r, er_r_retino, "BU_R") #Retinotopical
tb_er_r_retino.make_type_plots(plot_name = "Type Connections/Retinotopic TuBu_R to ER_R")


tb_er_l_retino = mapping.ConnectionMap(tb_l, er_l_retino, "BU_L") #Retinotopical
tb_er_l_retino.make_type_plots(plot_name = "Type Connections/Retinotopic TuBu_L to ER_L")






-------------
Compare Plots
-------------
neuread.ratio_plot(neuread.Ratio.TUBU_TO_METU, plot_name=\
    "Neuprint Comparison/Ratio of TuBu Neurons to MeTu Neurons per Hemisphere",\
    y_ticks=np.arange(0,0.8,0.2))

    
neuread.ratio_plot(neuread.Ratio.RING_TO_TUBU, plot_name=\
    "Neuprint Comparison/Ratio of Ring Neurons to TuBu Neurons per Hemisphere",\
    y_ticks=np.arange(0,4,1))

    
neuread.lobula_counts(just_metu4=False,\
    plot_name="Neuprint Comparison/Lobula Synapse Count in Hemibrain MeTu")


-------------
Miscellaneous
-------------
specific.fafb_syn_counts(["MeTu4a_L", "MeTu4a_R", "MeTu4b_L", "MeTu4b_R", 
    "MeTu4c_L", "MeTu4c_R", "MeTu4d_L", "MeTu4d_R"], "LO", \
    plot_name="Miscellaneous/Lobula Synapse Count per MeTu4 Neuron",
    y_axis="Lobula Synapse Count")
    
    
specific.tutu_comparison()





---------------
Synapse Density
---------------
colors = ["TPurple", "TPurple", "TBlue", "TRed", "TGreen", "TOrange2", 
          "TOrange2", "TRed", "TGreen", "TPink"]
for i, j in zip(maps, colors):
    density.plot_maps(maps[i], color=j, plot_name=f"{i} Blurrier", blur=10, 
                      max_value=0.157)

    
"""
"""
-----------
Maintenance
-----------


maintenance.find_undone_input()
================
#Last Run 5/8/23, identified 9 neurons. 


maintenance.check_joined_concurrence(all_mt_l + all_mt_r + tb_l + tb_r + \
                                     er_all + aotu46 + tutu + upstream_neurons + \
                                     ["Mi1_L", "Mi1_R"])
================
#Last Run 8/23/23
No issues

maintenance.check_proof_concurrence()
================
#Last Run 6/6/23
No issues


maintenance.check_duplicates()


To find ID discrepancies:
mapping.add_types(mt_l + mt_r)
mapping.add_types(tb_l + tb_r + ["TuBu_misc_R"])
"""
"""
Colors
======
"#FF0000" - Red for presynaptic
"#00FFFF" - Cyan for postsynaptic

"#D62728" - Red for compare figure and MeTu1
"#1F77B4" - Blue for compare figure and MeTu2
"#2CA02C" - Green for compare figure and MeTu3
"#FF8F0E" - Better Yellow for compare figure and MeTu4
"#FF7F0E" - Orange
"#BC27D6" - Purple
"#D627B5" - Pink
"#27B0D6" - Cyan
"#FFA90E" - Old Yellow for compare figure

Outdated Colors
===============
"#A0451C" - Dark Red for PL neurons
"#1F3481" - Dark Blue for PC neurons
"#477227" - Dark Green for A neurons
"#BB8627" - Dark Yellow for M neurons


"""
"""
Making URLs
-----------
pl_l = mapping.neur_ids["TuBu08_L"]
pc_l = np.concatenate((mapping.neur_ids["TuBu01_L"], mapping.neur_ids["TuBu06_L"]))
a_l = np.concatenate((mapping.neur_ids["TuBu07_L"], mapping.neur_ids["TuBu09_L"], \
                      mapping.neur_ids["TuBu10_L"]))
m_l = np.concatenate((mapping.neur_ids["TuBu02_L"], mapping.neur_ids["TuBu03_L"], \
                      mapping.neur_ids["TuBu04_L"], mapping.neur_ids["TuBu05_L"]))

neurs = np.concatenate((pl_l, pc_l, a_l, m_l))
colors = ["#A0451C"] * len(pl_l) + ["#1F3481"] * len(pc_l) + ["#477227"] * len(a_l) \
        + ["#BB8627"] * len(m_l)

flywire.encode_url(segments=neurs, open_browser=True, seg_colors=colors)





url_mt_l = np.concatenate(([mapping.neur_ids[x] for x in mt_l]))
url_mt_r = np.concatenate(([mapping.neur_ids[x] for x in mt_r]))
url_tb_l = np.concatenate(([mapping.neur_ids[x] for x in tb_l]))
url_tb_r = np.concatenate(([mapping.neur_ids[x] for x in tb_r]))
url_er_l = np.concatenate(([mapping.neur_ids[x] for x in er_l]))
url_er_r = np.concatenate(([mapping.neur_ids[x] for x in er_r]))

neurs = np.concatenate((url_mt_l, url_mt_r, url_tb_l, url_tb_r, url_er_l, url_er_r))
colors = ["#D627B5"] * len(url_mt_l) + ["#BC27D6"] * len(url_mt_r) + \
            ["#1F77B4"] * len(url_tb_l) + ["#FF7F0E"] * len(url_tb_r) + \
            ["#D62728"] * len(url_er_l) + ["#2CA02C"] * len(url_er_r)
flywire.encode_url(segments=neurs, open_browser=True, seg_colors=colors)


"""
"""
Meshes
------

current_type = mt_r

mapping.add_types(current_type)
all_ids = np.array([], dtype = np.int64)
for i in current_type:
    all_ids = np.append(all_ids, mapping.neur_ids[i])
meshes.cv.mesh.save(all_ids)


"""


def find_path_old(output_region = "LAL_R", min_syn_count = 0):
    df = readfiles.import_file("Codex Connections", file_type = "csv")
    df = df[df["syn_count"] >= min_syn_count]
    aotu_df = df[df["neuropil"] == "AOTU_R"]
    lal_df = df[df["neuropil"] == output_region]
    post_aotu = np.unique(aotu_df["post_root_id"])
    pre_lal = np.unique(lal_df["pre_root_id"])
    
    aotu_neurs = np.array([], dtype = np.int64)
    for i in post_aotu:
        if not i in pre_lal:
            continue
        aotu_neurs = np.append(aotu_neurs, i)
    
    print("AOTU Neuron Count:", len(aotu_neurs))
    
    searched_neurs = np.array([], dtype = np.int64)
    OpticTu = np.array([], dtype = np.int64)
    for i in aotu_neurs:
        temp_df = aotu_df[aotu_df["post_root_id"] == i]
        pre_neurs = np.unique(temp_df["pre_root_id"])
        for j in pre_neurs:
            if j in searched_neurs:
                continue
            searched_neurs = np.append(searched_neurs, j)
            temp_df2 = df[df["post_root_id"] == j]
            temp_df2 = temp_df2[(temp_df2["neuropil"] == "ME_R") |
                                (temp_df2["neuropil"] == "LO_R") |
                                (temp_df2["neuropil"] == "LOP_R")]
            if len(temp_df2) == 0:
                continue
            OpticTu = np.append(OpticTu, j)
    
    print("Optic Lobe Neuron Count:", len(OpticTu))
    return OpticTu
            
def get_lc10():
    labels = readfiles.import_file("Codex Labels", file_type = "csv")
    lc10 = np.array([], dtype = np.int64)
    for i in ["Lobula columnar 10; Lc10", "LC10", "LC10 (Wu, 2016)"]:
        temp_df = labels[labels["label"] == i]
        lc10 = np.concatenate((lc10, np.asarray(temp_df["root_id"])))
    return np.unique(lc10)


def find_connections(output_region = "LAL_R", min_syn_count = 3):
    df = readfiles.import_file("Codex Connections", file_type = "csv")
    df = df[df["syn_count"] >= min_syn_count]
    aotu_df = df[df["neuropil"] == "AOTU_R"]
    lal_df = df[df["neuropil"] == output_region]
    post_aotu = np.unique(aotu_df["post_root_id"])
    pre_lal = np.unique(lal_df["pre_root_id"])
    
    aotu_neurs = np.array([], dtype = np.int64)
    for i in post_aotu:
        if not i in pre_lal:
            continue
        aotu_neurs = np.append(aotu_neurs, i)
    
    pre_neurs = np.array([], dtype = np.int64)
    #lc10 = get_lc10()
    for i in aotu_neurs:
        temp_df = aotu_df[(aotu_df["post_root_id"] == i) &
                          ~(aotu_df["pre_root_id"].isin(pre_neurs))
                          #& (aotu_df["pre_root_id"].isin(lc10))
                          ]
        for j in np.unique(temp_df["pre_root_id"]):
            temp_df2 = df[df["post_root_id"] == j]
            temp_df2 = temp_df2[(temp_df2["neuropil"] == "ME_R") |
                                (temp_df2["neuropil"] == "LO_R") |
                                (temp_df2["neuropil"] == "LOP_R")]
            if len(temp_df2) == 0:
                continue
            pre_neurs = np.append(pre_neurs, j)
    
    connections = np.zeros((len(pre_neurs), len(aotu_neurs)), dtype = int)
    for indi, i in enumerate(pre_neurs):
        for indj, j in enumerate(aotu_neurs):
            temp_df = aotu_df[(aotu_df["pre_root_id"] == i) &
                              (aotu_df["post_root_id"] == j)]
            temp_arr = np.asarray(temp_df["syn_count"])
            if len(temp_arr) == 0:
                continue
            connections[indi][indj] = temp_arr[0]
    
    labels = readfiles.import_file("Codex Labels", file_type = "csv")

    aotu_labels = np.array([], dtype=str)
    for i in aotu_neurs:
        temp_df = labels[labels["root_id"] == i]
        if len(temp_df) == 0:
            aotu_labels = np.append(aotu_labels, str(i))
            continue
        temp_arr = np.asarray(temp_df["label"], dtype=str)
        aotu_labels = np.append(aotu_labels, f"{i!s} temp_arr[0]")
    
    pre_labels = np.array([], dtype=str)
    for i in pre_neurs:
        temp_df = labels[labels["root_id"] == i]
        if len(temp_df) == 0:
            pre_labels = np.append(pre_labels, str(i))
            continue
        temp_arr = np.asarray(temp_df["label"], dtype=str)
        pre_labels = np.append(pre_labels, f"{i!s} temp_arr[0]")
        
    
    conn_df = pd.DataFrame(connections, columns = aotu_labels)
    #conn_df.insert(0, "", np.asarray(pre_neurs, dtype = str))
    conn_df.insert(0, "", pre_labels)
    
    return conn_df



def find_path_lc10_to_glno():
    mapping.add_types([f"{x}_{y}" for x in ["LC10","GLNO"] for y in "LR"])
    lc10_r = mapping.neur_ids["LC10_R"]
    glno = mapping.ids_from_types([f"GLNO_{x}" for x in "LR"])
    return find_path(lc10_r, glno)


def find_path(pre_neurs, goal_neurs, min_neurs=10, counter=0, 
              seen_neurs=np.array([])):
    print(f"Layer: {counter}, Number of Neurons: {len(pre_neurs)!s}")
    syn_df = fw.fetch_synapses(pre_neurs, post=False)
    temp_df = syn_df[syn_df["post"].isin(goal_neurs)]
    if len(temp_df!=0):
        return {counter: np.unique(temp_df["pre"]), 
                counter+1: np.unique(temp_df["post"])}
    next_layer = utils.count_instances(list(syn_df["post"]), min_neurs=min_neurs)
    seen_neurs = np.unique(np.concatenate((seen_neurs, np.array(pre_neurs))))
    for i in list(next_layer.keys()):
        if i in seen_neurs:
            del next_layer[i]
    path_end = find_path(np.unique(list(next_layer.keys())), goal_neurs, 
                         min_neurs=min_neurs,
                         counter=counter+1, seen_neurs=seen_neurs)
    next_layer_neurs = path_end[counter+1]
    syn_df = syn_df[syn_df["post"].isin(next_layer_neurs)]
    path_end[counter] = np.unique(syn_df["pre"])
    return path_end
    
    
"""
Most Prominent Path
-------------------
LC10_R
V
720575940632409500
V
720575940631043139 (183 connections to GLNO)
V
GLNO



"""    




    
"""
LC10_R
V
GLNO_Path_Layer_2
720575940604584416, 720575940604732350, 720575940605920480,
       720575940606610976, 720575940606739390, 720575940607125771,
       720575940607887746, 720575940608529710, 720575940608661955,
       720575940609234130, 720575940609808905, 720575940610389873,
       720575940610903961, 720575940611112818, 720575940611190852,
       720575940611633010, 720575940611675634, 720575940612031469,
       720575940612218547, 720575940612227880, 720575940612752611,
       720575940613189785, 720575940613502786, 720575940613588246,
       720575940613643122, 720575940613904638, 720575940614190022,
       720575940614311746, 720575940614543527, 720575940614671411,
       720575940615258178, 720575940615522922, 720575940616012061,
       720575940616102550, 720575940616169117, 720575940616169629,
       720575940616334386, 720575940616676418, 720575940617101713,
       720575940617151026, 720575940617260374, 720575940617520413,
       720575940617654341, 720575940617909643, 720575940618254750,
       720575940618380081, 720575940618578150, 720575940618628677,
       720575940618959517, 720575940619021272, 720575940619859547,
       720575940620176286, 720575940620239739, 720575940620352459,
       720575940620362049, 720575940620444788, 720575940620686964,
       720575940620836988, 720575940621195290, 720575940621371045,
       720575940621673852, 720575940621679275, 720575940621751174,
       720575940621779661, 720575940621988461, 720575940622029990,
       720575940622053185, 720575940622053229, 720575940622105182,
       720575940622143245, 720575940622364598, 720575940622400648,
       720575940622431860, 720575940622494699, 720575940622538520,
       720575940622680724, 720575940622706420, 720575940623004933,
       720575940623107509, 720575940623318280, 720575940623463126,
       720575940623545843, 720575940623688013, 720575940623827533,
       720575940624065279, 720575940624304423, 720575940624305407,
       720575940624310468, 720575940624382695, 720575940624575076,
       720575940624579100, 720575940624589252, 720575940624599492,
       720575940624607940, 720575940624856112, 720575940625019696,
       720575940625090823, 720575940625440100, 720575940625498376,
       720575940625651849, 720575940625980403, 720575940626050599,
       720575940626283066, 720575940626552676, 720575940626733374,
       720575940626960688, 720575940626990480, 720575940627068772,
       720575940627076964, 720575940627094696, 720575940627144138,
       720575940627231400, 720575940627240783, 720575940627279674,
       720575940627351140, 720575940627497244, 720575940627717246,
       720575940627777067, 720575940628119339, 720575940628120744,
       720575940628227198, 720575940628342282, 720575940628653255,
       720575940628873722, 720575940628916444, 720575940628936168,
       720575940628983784, 720575940629084524, 720575940629290191,
       720575940629532777, 720575940629539628, 720575940629883971,
       720575940630140119, 720575940630191447, 720575940630310431,
       720575940630528872, 720575940630788222, 720575940630826660,
       720575940630882263, 720575940630924791, 720575940630940943,
       720575940630954434, 720575940631253820, 720575940631263023,
       720575940631441372, 720575940631536415, 720575940631839457,
       720575940632102607, 720575940632409500, 720575940632871322,
       720575940632941881, 720575940632943277, 720575940633429331,
       720575940633573836, 720575940633601249, 720575940633616332,
       720575940633900896, 720575940634041756, 720575940634553230,
       720575940634743999, 720575940635154830, 720575940635385440,
       720575940635659872, 720575940635873946, 720575940636351845,
       720575940636491631, 720575940636589492, 720575940636759386,
       720575940637179759, 720575940637416078, 720575940637464718,
       720575940637846643, 720575940637901540, 720575940637952624,
       720575940638277488, 720575940638482421, 720575940638681203,
       720575940639031642, 720575940639318079, 720575940640175549,
       720575940640792539, 720575940641188304, 720575940641276045,
       720575940641780187, 720575940643302797, 720575940646407476,
       720575940650999417, 720575940652922529, 720575940653696502,
       720575940658655873
V
GLNO_Path_Layer_3
720575940606664894, 720575940607982345, 720575940608140380,
       720575940609109602, 720575940609284728, 720575940609405020,
       720575940610075346, 720575940612143075, 720575940613891442,
       720575940614441834, 720575940614461014, 720575940615101747,
       720575940615620002, 720575940616867382, 720575940617374242,
       720575940618722542, 720575940620027515, 720575940620431798,
       720575940621895006, 720575940623120222, 720575940623260201,
       720575940624287821, 720575940624347540, 720575940624394790,
       720575940625673964, 720575940625849236, 720575940627578055,
       720575940627831641, 720575940627849515, 720575940628034603,
       720575940628321352, 720575940628494695, 720575940628694019,
       720575940629415936, 720575940629683500, 720575940630881337,
       720575940630986999, 720575940631022071, 720575940631043139,
       720575940631690860, 720575940632126860, 720575940632380216,
       720575940633145903, 720575940633631396, 720575940637064926,
       720575940639278399, 720575940639278781, 720575940640608219,
       720575940640749939, 720575940640764733, 720575940653288097,
       720575940654643873
V
GLNO

"""
    










=======
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
import flycode.graphics as graphics
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



"""
--------
Figure 1
--------
#Fig. 1c
colors = ["TPurple", "TPurple", "TBlue", "TRed", "TGreen", "TOrange2", 
          "TOrange2", "TRed", "TGreen", "TPink"]
for i, j in zip(maps, colors):
    density.plot_maps(maps[i], color=j, plot_name=f"{i} Blurrier", blur=10, 
                      max_value=0.157)

#Fig. 1di
mt_tb_r = mapping.ConnectionMap(mt_r_general, tb_r, "AOTU_R")
mt_tb_r.make_type_plots(plot_name = "MeTu_R to TuBu_R",
                        fig_size = fig_size_mt_tb)

#Fig. 1dii
mt_tb_l = mapping.ConnectionMap(mt_l_general, tb_l, "AOTU_L")
mt_tb_l.make_type_plots(plot_name = "MeTu_L to TuBu_L",
                        fig_size = fig_size_mt_tb)

#Fig. 1diii
mt_tb_r = mapping.ConnectionMap(mt_r, tb_r, "AOTU_R") #Specific
mt_tb_r.make_type_plots(plot_name = "Subtyped MeTu_R to TuBu_R")

#Fig. 1div
mt_tb_l = mapping.ConnectionMap(mt_l, tb_l, "AOTU_L") #Specific
mt_tb_l.make_type_plots(plot_name = "Subtyped MeTu_L to TuBu_L")


---------
Figure S1
---------
#Fig. S1biii and S1biv
f1_scores = volumes.compare_f1_scores()
volumes.plot_f1(f1_scores)
volumes.get_all_t_tests(f1_scores) #Gets t-test p-scores.


---------
Figure S2
---------
#Figure this one out


---------
Figure S3
---------
#Fig. S3a
bihem_aotu_r = mapping.ConnectionMap(aotu46+tutu+mt_r, aotu46+tutu+mt_r+tb_r, "AOTU_R")
bihem_aotu_r.make_type_plots(\
            plot_name="Bihemispheric Connections (AOTU_R)",
            fig_size=bihem_size)

#Fig. S3b
bihem_aotu_l = mapping.ConnectionMap(aotu46+tutu+mt_l, aotu46+tutu+mt_l+tb_l, "AOTU_L")
bihem_aotu_l.make_type_plots(\
            plot_name="Bihemispheric Connections (AOTU_L)",
            fig_size=bihem_size)

#Fig. S3c
bihem_bu_r = mapping.ConnectionMap(aotu46+tb_r, aotu46+tb_r+er_r, "BU_R")
bihem_bu_r.make_type_plots(\
            plot_name="Bihemispheric Connections (BU_R)",
            fig_size=bihem_size)

#Fig. S3d
bihem_bu_l = mapping.ConnectionMap(aotu46+tb_l, aotu46+tb_l+er_l, "BU_L")
bihem_bu_l.make_type_plots(\
            plot_name="Bihemispheric Connections (BU_L)",
            fig_size=bihem_size)


---------
Figure S4
---------
#Fig. S4
eb_inter = mapping.ConnectionMap(er_all_plus_lal+exr+eb_broad, 
                                 er_all_plus_lal+exr+eb_broad, "EB")
eb_inter.make_type_plots(plot_name = "EB Interconnectivity",
                         fig_size = (8, 8))


---------
Figure S5
---------
#Fig. S5ai
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_me = mapping.ConnectionMap(mt1_r, mt1_r, "ME_R", exclude=mt1_r_outliers)
mt1_r_inter_me.make_connectivity_plots(plot_name="MT1_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt1)

#Fig. S5aii
mt1_l_inter_me = mapping.ConnectionMap(mt1_l, mt1_l, "ME_L")
mt1_l_inter_me.make_connectivity_plots(plot_name="MT1_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt1)

#Fig. S5bi
mt1_r_outliers = fw.locs_to_segments(mt1_r_outlier_coords)
mt1_r_inter_aotu = mapping.ConnectionMap(mt1_r, mt1_r, "AOTU_R",)
mt1_r_inter_aotu.make_connectivity_plots(plot_name="MT1_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt1)

#Fig. S5bii
mt1_l_inter_aotu = mapping.ConnectionMap(mt1_l, mt1_l, "AOTU_L")
mt1_l_inter_aotu.make_connectivity_plots(plot_name="MT1_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt1)

#Fig. S5ci
mt_tb_pl_r = mapping.ConnectionMap(mt1_r, tb_pl_r, "AOTU_R")
mt_tb_pl_r.make_connectivity_plots(plot_name="MeTu1_R to TuBu_R",
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S5cii
mt_tb_pl_l = mapping.ConnectionMap(mt1_l, tb_pl_l, "AOTU_L")
mt_tb_pl_l.make_connectivity_plots(plot_name="MeTu1_L to TuBu_L",
                                   fig_size=fig_size_mt_tb_pl)

#Fig. S5di
tb_er_pl_r = mapping.ConnectionMap(tb_pl_r, er_pl_r, "BU_R")
tb_er_pl_r.make_connectivity_plots(plot_name="TuBu_PL_R to ER_PL_R",
                                   fig_size=fig_size_tb_er_pl_r)

#Fig. S5dii
fig_size_tb_er_pl_l = (1.15, 1.15)
tb_er_pl_l = mapping.ConnectionMap(tb_pl_l, er_pl_l, "BU_L")
tb_er_pl_l.make_connectivity_plots(plot_name="TuBu_PL_L to ER_PL_L",
                                   fig_size=fig_size_tb_er_pl_l)


---------
Figure S6
---------
#Fig. S6b, S6di-ii
specific.make_bihem_pie_charts()

#Fig. S6c
comparison.nt_by_types(["AOTU046_L_1", "AOTU046_L_2", "AOTU046_R_1", "AOTU046_R_2"],
        palette=["#e6194B", "#EBC400", "#3E7748", "#4400dd"],
        plot_names=["AOTU046 NT Predictions"], fig_size=(1, 0.75))

#Fig. S6ei
comparison.nt_by_types(["TuTuB_a_L", "TuTuB_a_R"], 
                       plot_names=["TuTuB_a NT Predictions"],
                       palette=["#e6194B", "#3E7748"])

#Fig. S6eii
comparison.nt_by_types(["TuTuB_b_L", "TuTuB_b_R"],
                       plot_names=["TuTuB_b NT Predictions"],
                       palette=["#e6194B", "#3E7748"])

#Fig. S3f
specific.tutu_comparison(plot_name="TuTu Synapse Counts Comparison")


---------
Figure S7
---------
#Fig. S7ai
mt2_r_outliers = fw.locs_to_segments(mt2_r_outlier_coords)
mt2_r_inter_me = mapping.ConnectionMap(mt2_r, mt2_r, "ME_R", exclude = mt2_r_outliers)
mt2_r_inter_me.make_connectivity_plots(plot_name="MT2_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt2)

#Fig. S7aii
mt2_l_inter_me = mapping.ConnectionMap(mt2_l, mt2_l, "ME_L")
mt2_l_inter_me.make_connectivity_plots(plot_name="MT2_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt2)

#Fig. S7bi
mt2_r_inter_aotu = mapping.ConnectionMap(mt2_r, mt2_r, "AOTU_R",)
mt2_r_inter_aotu.make_connectivity_plots(plot_name="MT2_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt2)

#Fig. S7bii
mt2_l_inter_aotu = mapping.ConnectionMap(mt2_l, mt2_l, "AOTU_L")
mt2_l_inter_aotu.make_connectivity_plots(plot_name="MT2_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt2)

#Fig. S7ci
mt_tb_pc_r = mapping.ConnectionMap(mt2_r, tb_pc_r, "AOTU_R")
mt_tb_pc_r.make_connectivity_plots(plot_name="MeTu2_R to TuBu_R",
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S7cii
mt_tb_pc_l = mapping.ConnectionMap(mt2_l, tb_pc_l, "AOTU_L")
mt_tb_pc_l.make_connectivity_plots(plot_name="MeTu2_L to TuBu_L",
                                   fig_size=fig_size_mt_tb_pc)

#Fig. S7di
tb_er_pc_r = mapping.ConnectionMap(tb_pc_r, er_pc_r, "BU_R")
tb_er_pc_r.make_connectivity_plots(plot_name="TuBu_PC_R to ER_PC_R",
                                   fig_size=fig_size_tb_er_pc)

#Fig. S7dii
tb_er_pc_l = mapping.ConnectionMap(tb_pc_l, er_pc_l, "BU_L")
tb_er_pc_l.make_connectivity_plots(plot_name="TuBu_PC_L to ER_PC_L",
                                   fig_size=fig_size_tb_er_pc)


---------
Figure S8
---------
#Fig. S8ai
mt3_r_inter_me = mapping.ConnectionMap(mt3_r, mt3_r, "ME_R")
mt3_r_inter_me.make_connectivity_plots(plot_name="MT3_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt3)

#Fig. S8aii
mt3_l_inter_me = mapping.ConnectionMap(mt3_l, mt3_l, "ME_L")
mt3_l_inter_me.make_connectivity_plots(plot_name="MT3_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt3)

#Fig. S8bi
mt3_r_inter_aotu = mapping.ConnectionMap(mt3_r, mt3_r, "AOTU_R")
mt3_r_inter_aotu.make_connectivity_plots(plot_name="MT3_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt3)

#Fig. S8bii
mt3_l_inter_aotu = mapping.ConnectionMap(mt3_l, mt3_l, "AOTU_L")
mt3_l_inter_aotu.make_connectivity_plots(plot_name="MT3_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt3)

#Fig. S8ci
mt_tb_a_r = mapping.ConnectionMap(mt3_r, tb_a_r, "AOTU_R")
mt_tb_a_r.make_connectivity_plots(plot_name="MeTu3_R to TuBu_R",
                                  fig_size=fig_size_mt_tb_a)

#Fig. S8cii
mt_tb_a_l = mapping.ConnectionMap(mt3_l, tb_a_l, "AOTU_L")
mt_tb_a_l.make_connectivity_plots(plot_name="MeTu3_L to TuBu_L",
                                  fig_size=fig_size_mt_tb_a)

#Fig. S8di
tb_er_a_r = mapping.ConnectionMap(tb_a_r, er_a_r, "BU_R")
tb_er_a_r.make_connectivity_plots(plot_name="TuBu_A_R to ER_A_R",
                                  fig_size=fig_size_tb_er_a)

#Fig. S8dii
tb_er_a_l = mapping.ConnectionMap(tb_a_l, er_a_l, "BU_L")
tb_er_a_l.make_connectivity_plots(plot_name="TuBu_A_L to ER_A_L",
                                  fig_size=fig_size_tb_er_a)

#Fig. S8e
specific.sm17_map_by_dv_axis(plot_name= 
        "MeTu3bc Vertical Positioning by Connections to Sm17")

#Fig. S8g
specific.mt3_pre_connections(["Sm23", "MeMeDRA", "Mi15"], 
                             plot_name="MeTu3 Presynaptic")


---------
Figure S9
---------
#Fig. S9ai
mt4_r_outliers = fw.locs_to_segments(mt4_r_outlier_coords)
mt4_r_inter_me = mapping.ConnectionMap(mt4_r, mt4_r, "ME_R", exclude = mt4_r_outliers)
mt4_r_inter_me.make_connectivity_plots(plot_name="MT4_R Interconnectivity (ME_R)",
                                       fig_size=fig_size_mt4)

#Fig. S9aii
mt4_l_outliers = fw.locs_to_segments(mt4_l_outlier_coords)
mt4_l_inter_me = mapping.ConnectionMap(mt4_l, mt4_l, "ME_L",)
mt4_l_inter_me.make_connectivity_plots(plot_name="MT4_L Interconnectivity (ME_L)",
                                       fig_size=fig_size_mt4)

#Fig. S9bi
mt4_r_inter_aotu = mapping.ConnectionMap(mt4_r, mt4_r, "AOTU_R",)
mt4_r_inter_aotu.make_connectivity_plots(plot_name="MT4_R Interconnectivity (AOTU_R)",
                                         fig_size=fig_size_mt4)

#Fig. S9bii
mt4_l_inter_aotu = mapping.ConnectionMap(mt4_l, mt4_l, "AOTU_L",)
mt4_l_inter_aotu.make_connectivity_plots(plot_name="MT4_L Interconnectivity (AOTU_L)",
                                         fig_size=fig_size_mt4)

#Fig. S9ci
mt_tb_m_r = mapping.ConnectionMap(mt4_r, tb_m_r, "AOTU_R")
mt_tb_m_r.make_connectivity_plots(plot_name="MeTu4_R to TuBu_R",
                                  fig_size=fig_size_mt_tb_m)

#Fig. S9cii
mt_tb_m_l = mapping.ConnectionMap(mt4_l, tb_m_l, "AOTU_L")
mt_tb_m_l.make_connectivity_plots(plot_name="MeTu4_L to TuBu_L",
                                  fig_size=fig_size_mt_tb_m)

#Fig. S9di
tb_er_m_r = mapping.ConnectionMap(tb_m_r, er_m_r, "BU_R")
tb_er_m_r.make_connectivity_plots(plot_name="TuBu_M_R to ER_M_R",
                                  fig_size=fig_size_tb_er_m)

#Fig. S9dii
tb_m_l.insert(1, "TuBu_misc_L")
tb_er_m_l = mapping.ConnectionMap(tb_m_l, er_m_l, "BU_L")
tb_er_m_l.make_connectivity_plots(plot_name="TuBu_M_L to ER_M_L",
                                  fig_size=fig_size_tb_er_m)
tb_m_l.remove("TuBu_misc_L")


----------
Figure S12
----------
#Fig. S12a
neuread.plot_comparison(gen_metu,
        datasets={"FAFB": "LR", "Hemibrain": "R", "FIB-SEM": "R"},
        plot_name="MeTu Neuron Counts Comparison", 
        y_ticks=np.arange(0, 140, 20))

#Fig. S12b
neuread.plot_comparison(gen_tubu, 
        plot_name="TuBu Neuron Counts Comparison", 
        y_ticks=np.arange(0, 14, 2))

#Fig. S12c
neuread.plot_comparison(gen_ring,
        datasets={"FAFB": "LR", "Hemibrain": "LR"},
        plot_name="Ring Neuron Counts Comparison", 
        y_ticks=np.arange(0, 14, 2), fig_size = (2.4, 1.25))

#Fig. S12d
neuread.plot_comparison(gen_bihem, 
        datasets={"FAFB": "LR", "Hemibrain": "LR"},
        plot_name="Bihemispheric Neuron Counts Comparison", 
        y_ticks=np.arange(0, 3, 1))

#Fig. S12e
neuread.ratio_plot(neuread.Ratio.METU_TO_TUBU, plot_name=\
    "Ratio of MeTu Neurons to TuBu Neurons per Hemisphere",\
    y_ticks=np.arange(0, 18, 3))

#Fig. S12f
neuread.ratio_plot(neuread.Ratio.TUBU_TO_RING, plot_name=\
    "Ratio of TuBu Neurons to Ring Neurons per Hemisphere",\
    y_ticks=np.arange(0,5,1))

#Fig. S12g
neuread.lobula_counts(just_metu4=True,\
        plot_name="Lobula Synapse Count in Hemibrain MeTu4")

#Fig. S12h
specific.full_comparison([f"MeTu4{x}" for x in "abcd"], "LO", 
        plot_name="Lobula Synapse Count per MeTu4 Neuron", 
        y_axis="Lobula Synapse Count")

#Fig. S12i_i-ii
comparison.mt4e_partner_comparison("Fig S1/MeTu4e/Partner Type Comparison", True)

#Fig. S12i_iii-viii
comparison.mt4e_dorsal_comparison("Fig S1/MeTu4e/Comparing", True)


----------
Figure S13
----------
#Fig. S13ai
comparison.compare_metu("Medulla Presynapse Counts")

#Fig. S13aii
comparison.compare_metu("Medulla Postsynapse Counts")

#Fig. S13aiii
comparison.compare_metu("Column Counts")

#Fig. S13aiv
comparison.compare_metu("AOTU Presynapse Counts")

#Fig. S13av
comparison.compare_metu("AOTU Postsynapse Counts")

#Fig. S13avi
comparison.compare_metu("Ellipse Ratio")

#Fig. S13bi-viii
for i in ["Medulla Postsynapse Count", "AOTU Postsynapse Count",
          "Medulla Presynapse Count", "AOTU Presynapse Count",
          "Column Count", "Ellipse Ratio", "Ellipse Major Axis Length (nm)",
          "Ellipse Minor Axis Length (nm)"]:
    for j in ["A-P", "V-D"]:
        comparison.scatter_plots(x=f"Relative Offset from Medulla Centroid ({j})",
                                 y=i, plot_folder=f"{i} {j}/")


----------
Figure S14
----------
#Fig. S14a-j
comparison.nt_by_types(
    [f"{x}_R" for x in comparison.colors], 
    palette=list(comparison.colors.values()), 
    plot_names=[f"Fig S6/Neurotransmitters/{x}" for x in comparison.colors],
    save_figure=True,
    fig_size=(1,0.75), 
    separate_plots=True)


----------
Figure S15
----------
#Fig. S15c
ring_to_epg = mapping.ConnectionMap(er_all, epg, "EB")
ring_to_epg.make_type_plots("Fig S30/Ring to EPG", fig_size=(3.0, 4.0))


------------
Spreadsheets
------------

#Download "Labels" from https://codex.flywire.ai/api/download and store it
#in "flycode/Readable" as "Codex Labels 783.csv"
label_df = proofreading.all_codex_names(all_in_paper)
utils.write_excel(label_df, "Codex Naming History")

percent_df = proofreading.full_percent_spread(all_in_paper)
utils.write_excel(percent_df, "FlyWire Consortium Edit Record")


------
Colors
------
"#FF0000" - Red for presynaptic
"#00FFFF" - Cyan for postsynaptic

"#D62728" - Red for compare figure and MeTu1
"#1F77B4" - Blue for compare figure and MeTu2
"#2CA02C" - Green for compare figure and MeTu3
"#FF8F0E" - Better Yellow for compare figure and MeTu4
"#FF7F0E" - Orange
"#BC27D6" - Purple
"#D627B5" - Pink
"#27B0D6" - Cyan
"""










>>>>>>> Stashed changes
