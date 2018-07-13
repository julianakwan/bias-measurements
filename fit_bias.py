#Script to make scatter plots and fit polynomial to delta_h = b_1*delta_m + b_2*delta_m^2 + eps
#Run it as python scatter_plot.py bin_no

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import sys
from scipy.optimize import curve_fit


def func(x,b1,b2,eps):
    return b1*x+b2*x*x+eps


if __name__ == "__main__":
    ngrid = 4
    step_no = 499
    scale_factor = 1.0
    mass_cut = np.arange(13.0,15.0,0.1)
    
    for model_no in range (19,20):
        nbins = 20
        b1 = np.zeros(nbins,dtype='float64')
        b2 = np.zeros(nbins,dtype='float64')
        eps = np.zeros(nbins,dtype='float64')

        #Use this for L2100 boxes
#        dens_m_filename = "/Users/jkwan/emulators/MiraU/density/M{:0>3d}/STEP{:d}/dens_M{:0>3d}_mm_L2100_N128_R4_{:5.4f}.dat".format(model_no, step_no, model_no, scale_factor)

        #Use this for L5000 boxes
        dens_m_filename = '/Users/jkwan/emulators/MiraU/density/M{:0>3d}/STEP{:d}dens_M{:0>3d}_mm_L5000_N256_R3.5_{:5.4f}.dat'.format(model_no, step_no, model_no, scale_factor)
        dens_m = np.fromfile(dens_m_filename, dtype='>f8')
        for binno in range (0,nbins):

            dens_h_filename = "/Users/jkwan/emulators/MiraU/density/M{:0>3d}/STEP{:d}/dens_M{:0>3d}_hh_L5000_N256_R3.5_sod_{:5.4f}_{:4.2f}.dat".format(model_no, step_no, model_no, scale_factor, mass_cut[binno])
            dens_h = np.fromfile(dens_h_filename, dtype='>f8')


            popt,pcov = curve_fit(func,dens_m,dens_h)
            b1[binno] = popt[0]  
            b2[binno] = popt[1]
            eps[binno] = popt[2]

            print(b1[binno])

        output_fname = "/Users/jkwan/emulators/bias/bias_measurements/STEP{:d}/M{:0>3d}_{:d}_sod_L5000_test.dat".format(step_no, model_no, step_no)
        np.savetxt(output_fname,np.transpose((mass_cut,b1,b2,eps)),fmt='%.4f')






