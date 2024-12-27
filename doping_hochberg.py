"""
This file creates a custom doping distribution for the 
"""
import tidy3d as td


def generate_doping_boxes():
    acceptor_boxes = []
    donors_boxes = []

    acceptor_boxes.append(td.ConstantDoping(box_coords=[[-5, 0, 0],[5, 0.22, 0]], concentration=1e15))

    # p implant
    acceptor_boxes.append(td.GaussianDoping(
        box_coords=[[-6, -0.3, 0],[-0.15, 0.098, 0]], 
        concentration=7e17,
        ref_con=1e6,
        width=0.1,
        source="ymax"))
    
    # n implant
    donors_boxes.append(td.GaussianDoping(
        box_coords=[[0.15, -0.3, 0],[6, 0.098, 0]], 
        concentration=5e17,
        ref_con=1e6,
        width=0.1,
        source="ymax") )
    
    # p++
    acceptor_boxes.append(td.GaussianDoping(
        box_coords=[[-6, -0.3, 0],[-2, 0.22, 0]], 
        concentration=1e19,
        ref_con=1e6,
        width=0.1,
        source="ymax"))
    
    # n++
    donors_boxes.append(td.GaussianDoping(
        box_coords=[[2, -0.3, 0],[6, 0.22, 0]], 
        concentration=1e19,
        ref_con=1e6,
        width=0.1,
        source="ymax"))
    
    
    # p wg implant
    acceptor_boxes.append(td.GaussianDoping(
        box_coords=[[-0.3, 0, 0],[0.06, 0.255, 0]], 
        concentration=5e17,
        ref_con=1e6,
        width=0.12,
        source="xmin"))
    
    # n wg implant
    donors_boxes.append(td.GaussianDoping(
        box_coords=[[-0.06, 0.02, 0],[0.25, 0.26, 0]], 
        concentration=7e17,
        ref_con=1e6,
        width=0.11,
        source="xmax"))
    
    return acceptor_boxes, donors_boxes