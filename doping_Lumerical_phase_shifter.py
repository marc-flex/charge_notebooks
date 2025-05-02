import numpy as np
from scipy.special import gamma, loggamma
import tidy3d as td

def pearson_type_iv_pdf(x, mu, beta, m, nu):
    """
    Compute the PDF of the Pearson Type IV distribution.
    """
    # Precompute constants
    const = np.abs(gamma(m + (nu / 2) * 1j))**2 / (gamma(m) * np.sqrt(np.pi * beta))
    term1 = (1 + (x - mu)**2 / beta)**(-m)
    term2 = np.exp(-nu * np.arctan((x - mu) / np.sqrt(beta)))
    
    return const * term1 * term2

def apply_pearson4(X: np.array, Y: np.array, doping_values: np.array, x_lim: list, y_lim: list, source_face: list, con: float, pearson_params: list):
    """
    X & Y are the coordinates of the SpatialDataArray to be filled out
    doping_values: the data of the SpatialDataArray
    x_lim: [x_min, x_max]
    y_lim: [y_min, y_max]
    source_face: list of zeros and ones specifying the source face. The source face is specified with a 1.
                The order is as follows: [x_min, x_max, y_min, y_max]
    pearson_params: list of parameters as they appear on Lumerical: range, straggle, skewness, kurtosis, lateral scatter
    NOTE: length units need to be consistent
    """
    # Pearson params
    rang = pearson_params[0]
    straggle = pearson_params[1]
    skewness = pearson_params[2]
    kurt = pearson_params[3]
    s = pearson_params[4] # lateral scatter
    
    mu = rang    # Location parameter - same as range in Lumerical
    m = 6/(kurt-3) + 2.5    # Shape parameter - approx 6/(kurtosis-3) + 2.5
    beta = straggle**2 * (2*m-3)  # Scale parameter - same is straggle**2*(2m-3)
    nu = skewness*np.sqrt(m)    # Skewness parameter

    # get max
    eta = np.linspace(-10, 10, 1000)
    valsp4 = pearson_type_iv_pdf(eta, mu, beta, m, nu)
    max_val = np.max(valsp4)

    if source_face[0] == 1 or source_face[1] == 1:
        for i, x in enumerate(X):
            x_contrib = 1
            if x >= x_lim[0] and x <= x_lim[1]:
                if source_face[0] == 1:
                    x_contrib = pearson_type_iv_pdf(x - x_lim[0], mu, beta, m, nu)
                if source_face[1] == 1:
                    x_contrib = pearson_type_iv_pdf(x_lim[1] - x, mu, beta, m, nu)
            else:
                x_contrib = 0
            
            for j, y in enumerate(Y):
                y_contrib = 1
                if y <= y_lim[0]:
                    y_contrib = np.exp(-(y - y_lim[0])*(y - y_lim[0])/s/s)
                if y >= y_lim[1]:
                    y_contrib = np.exp(-(y - y_lim[1])*(y - y_lim[1])/s/s)
                
                doping_values[i][j][0] = doping_values[i][j][0] + con * x_contrib * y_contrib / max_val

    if source_face[2] == 1 or source_face[3] == 1:
        for j, y in enumerate(Y):
            y_contrib = 1
            if y >= y_lim[0] and y <= y_lim[1]:
                if source_face[2] == 1:
                    y_contrib = pearson_type_iv_pdf(y - y_lim[0], mu, beta, m, nu)
                if source_face[3] == 1:
                    y_contrib = pearson_type_iv_pdf(y_lim[1] - y, mu, beta, m, nu)
            else:
                y_contrib = 0
            
            for i, x in enumerate(X):
                x_contrib = 1
                if x <= x_lim[0]:
                    x_contrib = np.exp(-(x - x_lim[0])*(x - x_lim[0])/s/s)
                if x >= x_lim[1]:
                    x_contrib = np.exp(-(x - x_lim[1])*(x - x_lim[1])/s/s)
                    
                doping_values[i][j][0] = doping_values[i][j][0] + con * x_contrib * y_contrib / max_val
    
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


def apply_constant(X: np.array, Y: np.array, doping_values: np.array, x_lim: list, y_lim: list, con: float):
    for i, x in enumerate(X):
        x_contrib = 0
        if x >= x_lim[0] and x <= x_lim[1]:
            x_contrib = 1
        
        for j, y in enumerate(Y):
            y_contrib = 0
            if y >= y_lim[0] and y <= y_lim[1]:
                y_contrib = 1

            doping_values[i][j][0] = doping_values[i][j][0] + con * x_contrib * y_contrib
    


Nx = 1000
Ny = 100
Nz = 1

delta = 0.1
X = np.linspace(-1.5 -delta, 1.55+delta, Nx)
Y = np.linspace(-0.05-delta, 0.4+delta, Ny)
acceptors = np.ones((Nx, Ny, Nz)) * 1e15
donors = np.zeros((Nx, Ny, Nz))

# wg p+
pearson_params = [0.3, 0.5, 0.5, 7, 0.1]
con = 2.5e17
source = [0, 0, 0, 1]
x_lim = [0.05, 1]
y_lim = [-0.5, 0.4]
apply_pearson4(X, Y, acceptors, x_lim, y_lim, source, con, pearson_params)

# wg n
pearson_params = [0.3, 0.5, 0.5, 7, 0.1]
con = 1.5e17
source = [0, 0, 0, 1]
x_lim = [-1, -0.05]
y_lim = [-0.5, 0.4]
apply_pearson4(X, Y, donors, x_lim, y_lim, source, con, pearson_params)

# source nwell
x_lim = [-1.75, -0.55]
y_lim = [-0.3, 0.001]
source_face = [0, 0, 0, 1]
con = 1e18
ref_con = 1e15
width = 0.12
apply_gaussian(X, Y, donors, x_lim, y_lim, source_face, con, ref_con, width)

# drain pwell
x_lim = [0.55, 1.75]
y_lim = [-0.3, 0.001]
source_face = [0, 0, 0, 1]
con = 1e18
ref_con = 1e14
width = 0.12
apply_gaussian(X, Y, acceptors, x_lim, y_lim, source_face, con, ref_con, width)

a = td.SpatialDataArray(data=acceptors-donors, coords={"x":X, "y":Y, "z": [0]})
d = td.SpatialDataArray(data=donors, coords={"x":X, "y":Y, "z": [0]})

import matplotlib.pyplot as plt
a.plot(y = "y")
plt.show()