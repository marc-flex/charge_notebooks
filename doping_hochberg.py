"""
This file creates a custom doping distribution for the 
"""
from tidy3d import SpatialDataArray
import numpy as np
import matplotlib.pyplot as plt


def apply_gaussian(X: np.array, Y: np.array, doping_values: np.array, x_lim: list, y_lim: list, source_face: list, con: float, ref_con: float, width: float):
    """
    X & Y are the coordinates of the SpatialDataArray to be filled out
    doping_values: the data of the SpatialDataArray
    x_lim: [x_min, x_max]
    y_lim: [y_min, y_max]
    source_face: list of zeros and ones specifying the source face. The source face is specified with a 1.
                The order is as follows: [x_min, x_max, y_min, y_max]
    NOTE: length units need to be consistent
    """
    # compute sigma value for the gaussian
    s = np.sqrt(-width*width/2/np.log(ref_con/con))
    for i, x in enumerate(X):
        x_contrib = 1
        if x >= x_lim[0] and x <= x_lim[0] + width and source_face[0] == 0:
            x_contrib = np.exp(-(x-x_lim[0]-width)*(x-x_lim[0]-width)/2/s/s)
        if x <= x_lim[1] and x >= x_lim[1]  -width and source_face[1] == 0:
            x_contrib = np.exp(-(x-x_lim[1]+width)*(x-x_lim[1]+width)/2/s/s)
        if x < x_lim[0] or x > x_lim[1]:
            x_contrib = 0
        
        for j, y in enumerate(Y):
            y_contrib = 1
            if y >= y_lim[0] and y <= y_lim[0] + width and source_face[2] == 0:
                y_contrib = np.exp(-(y-y_lim[0]-width)*(y-y_lim[0]-width)/2/s/s)
            if y <= y_lim[1] and y >= y_lim[1]  -width and source_face[3] == 0:
                y_contrib = np.exp(-(y-y_lim[1]+width)*(y-y_lim[1]+width)/2/s/s)
            if y < y_lim[0] or y > y_lim[1]:
                y_contrib = 0

            doping_values[i][j][0] = doping_values[i][j][0] + con * x_contrib * y_contrib

def get_dopings_hochberg():
    # all units in um
    # create a cartesian 2D mesh
    Nx = 2000
    Ny = 100
    Nz = 1

    delta = 0.1
    X = np.linspace(-5 -delta, 5+delta, Nx)
    Y = np.linspace(0-delta, 0.22+delta, Ny)
    acceptors = np.ones((Nx, Ny, Nz)) * 1e15
    donors = np.zeros((Nx, Ny, Nz))
        
    
    # p implant
    x_lim = list(np.array([-6, -0.15]))
    y_lim = list(np.array([-0.3, 0.098]))
    source_face = [0, 0, 0, 1]
    con = 7e17
    ref_con = 1e6
    width = 0.1
    apply_gaussian(X,Y,acceptors, x_lim, y_lim, source_face, con, ref_con, width)
    
    # n implant
    x_lim = list(np.array([0.15, 6]))
    y_lim = list(np.array([-0.3, 0.098]))
    source_face = [0, 0, 0, 1]
    con = 5e17
    ref_con = 1e6
    width = 0.1
    apply_gaussian(X,Y,donors, x_lim, y_lim, source_face, con, ref_con, width)
    
    # p++
    x_lim = list(np.array([-6, -2]))
    y_lim = list(np.array([-0.3, 0.22+delta]))
    source_face = [0, 0, 0, 1]
    con = 1e19
    ref_con = 1e6
    width = 0.1
    apply_gaussian(X,Y,acceptors, x_lim, y_lim, source_face, con, ref_con, width)
    
    # n++
    x_lim = list(np.array([2, 6]))
    y_lim = list(np.array([-0.3, 0.22+delta]))
    source_face = [0, 0, 0, 1]
    con = 1e19
    ref_con = 1e6
    width = 0.1
    apply_gaussian(X,Y,donors, x_lim, y_lim, source_face, con, ref_con, width)
    
    
    # p wg implant
    x_lim = list(np.array([-0.3, 0.06]))
    y_lim = list(np.array([0, 0.255]))
    source_face = [1, 0, 0, 0]
    con = 5e17
    ref_con = 1e6
    width = 0.12
    apply_gaussian(X,Y,acceptors, x_lim, y_lim, source_face, con, ref_con, width)
    
    # n wg implant
    x_lim = list(np.array([-0.06, 0.25]))
    y_lim = list(np.array([0.02, 0.26]))
    source_face = [0, 1, 0, 0]
    con = 7e17
    ref_con = 1e6
    width = 0.11
    apply_gaussian(X,Y,donors, x_lim, y_lim, source_face, con, ref_con, width)
    
    
    a = SpatialDataArray(data=acceptors, coords={"x":X, "y":Y, "z": [0]})
    d = SpatialDataArray(data=donors, coords={"x":X, "y":Y, "z": [0]})
    return a, d
    np.log(np.abs(a)).plot(y = "y")
    plt.show()