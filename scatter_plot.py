#Script to make scatter plots and fit polynomial to delta_h = b_1*delta_m + b_2*delta_m^2 + eps
#Run it as python scatter_plot.py bin_no

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import sys
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

def linear_bias(x,b1,eps):
    return b1*x+eps

def second_order_bias(x,b1,b2,eps):
    return b1*x+b2*x*x+eps


def third_order_bias(x,b1,b2,b3,eps):
    return b1*x+b2*x*x+b3**x*x*x+eps

def residuals(b,y,x):
    return y-(b[1]*x+b[2]*x*x+b[0])


def make_contour_plot(x, y, popt, resX=25, resY=25):
    xi = np.linspace(np.min(x), np.max(x), resX)
    yi = np.linspace(np.min(y), 2., resY)
    H, xi, yi = np.histogram2d(x,y,bins=(xi,yi))
    H = H.T
    fig = plt.figure()
    plt.imshow(H,interpolation='bicubic',origin='low',extent=[xi[0],xi[-1],yi[0],yi[-1]],aspect='auto',cmap='jet')
    plt.axis([xi[0],xi[-1],yi[0],yi[-1]])
    x = np.linspace(-1,np.max(x)+1,20)
    plt.plot(x,second_order_bias(x,*popt),'r-')
    plt.colorbar()
    plt.show()
    return fig


               

if __name__ == "__main__":
    binno = sys.argv[1]

    mass_cut = np.arange(12.5,15.0,0.1)

    print(mass_cut[int(binno)])

#For MiraU, dtype='>f8'. > = big endian, 8-bit float, i.e. double

    dens_m_filename = '/Users/jkwan/emulators/MiraU/density/M012/STEP499/dens_M012_mm_L5000_N256_R3.5_1.0000.dat'
    dens_m = np.fromfile(dens_m_filename, dtype='>f8')
    dens_h_filename = "/Users/jkwan/emulators/MiraU/density/M012/STEP499/dens_M012_hh_L5000_N256_R3.5_sod_1.0000_{0:4.2f}.dat".format(mass_cut[int(binno)])
    dens_h = np.fromfile(dens_h_filename, dtype='>f8')
     

#For TitanU:
#    dens_m = np.fromfile("/Users/jkwan/emulators/TitanU/density/M019/STEP247/dens_M019_mm_L2100_N128_R4_0.4985.dat", dtype='double')
#    dens_h_filename = "/Users/jkwan/emulators/TitanU/density/M019/STEP499/dens_M019_hh_L2100_N128_R4_sod_1.0000_{0:4.2f}.dat".format(mass_cut[int(binno)])
#    dens_h = np.fromfile(dens_h_filename, dtype='double')

    print(dens_h[0])
    print(dens_h[1])
    print(dens_m[0])
    
    popt,pcov = curve_fit(second_order_bias,dens_m,dens_h)


    print(popt[0])
    print(popt[1])
    print(popt[2])
    print(np.mean(dens_h))

    
    make_contour_plot(dens_m, dens_h, popt);
    
    index = np.random.choice(np.size(dens_m),10000)
    plt.scatter(dens_m[index],dens_h[index])


    x = np.linspace(-1,np.max(dens_m)+1,20)
    plt.plot(x,second_order_bias(x,*popt),'r-')
    plt.axis([-1,np.max(dens_m)+0.5,-1,np.max(dens_h)+0.5])
    plt.text(-0.8,0.9, r'b$_1$ = %.3f' % (popt[0]))
    plt.text(-0.8,0.8,  r'b$_2$ = %.3f' % (popt[1]))
    plt.xlabel(r"$\delta_m$", fontsize=14)
    plt.ylabel(r"$\delta_h$", fontsize=14)
    plt.show()

    



