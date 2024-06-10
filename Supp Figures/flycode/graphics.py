#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 14:22:24 2024

@author: dustin
"""

import numpy as np
import fpdf
import flycode.mapping as mapping


#0.040061633281972264
#0.04262877442273535
#0.08310502283105023
    

def lc10_to_glno_network():
    layers = [["LC10_R"], ["GLNO_Path_Layer_2"], ["GLNO_Path_Layer_3"], 
              ["GLNO_R", "GLNO_L"]]
    neurons = [mapping.ids_from_types(x) for x in layers]
    
    x_locations = np.linspace(20, 3180, 4, dtype=int)
    y_locations = [np.linspace(20, 7500, len(neurons[x]), dtype=int) for x in range(4)]
    
    pdf = fpdf.FPDF(orientation = 'P', unit = 'pt', format=(3200,7520))
    pdf.add_page()
    
    for i in range(3):
        pre = layers[i]
        post = layers[i+1]
        conn_map = mapping.ConnectionMap(pre, post, region="Connectome")
        weight_map = conn_map.weight_map
        print(np.max(weight_map))
        for indj, j in enumerate(weight_map):
            for indk, k in enumerate(j):
                if k==0:
                    continue
                width = (4/0.084)*k
                pdf.set_line_width(width)
                x1, x2 = x_locations[i], x_locations[i+1]
                y1, y2 = y_locations[i][indj], y_locations[i+1][indk]
                pdf.line(x1, y1, x2, y2)
    
    colors = [(214,39,40), (39,118,180), (44,160,44), (255,143,14)]
    for i in range(4):
        color = colors[i]
        pdf.set_fill_color(color[0],color[1],color[2])
        for j in y_locations[i]:
            pdf.ellipse(x_locations[i]-8, j-8, 16, 16, style="F")
    
    #Need to increment the figure name each time or it will generate the same figure
    pdf.output("/Users/dustin/Desktop/Test9.pdf", "F")
    pdf.close()
    
    return pdf




    
    
