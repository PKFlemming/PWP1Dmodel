import numpy as np

# of incident light i, i_sw is fraction presumed to be shortwave
# i_lw = 1 - i_sw is fr. presumed to be longwave

def fraction_reaching_depth_z(z, i_sw, lambda_lw, lambda_sw):
    i_lw = 1 - i_sw
    longwave_fraction_reaching_depth = i_lw * np.exp(-z / lambda_lw)
    shortwave_fraction_reaching_depth = i_sw * np.exp(-z / lambda_sw)
    fraction_reaching_depth = longwave_fraction_reaching_depth + shortwave_fraction_reaching_depth
    return fraction_reaching_depth

def fraction_absorbed(z, dz=1, i_sw = 0.6, lambda_lw = 20, lambda_sw = 0.62):
    fraction_reaching_top = fraction_reaching_depth_z(z, i_sw, lambda_lw, lambda_sw)
    fraction_reaching_bottom = fraction_reaching_depth_z(z+dz, i_sw, lambda_lw, lambda_sw)
    fraction_absorbed = fraction_reaching_top - fraction_reaching_bottom
    return fraction_absorbed
