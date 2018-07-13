import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import sys

mass = np.arange(13.0,15.0,0.1)


nbins = 20
nmodels = 26

bin_no = sys.argv[1]
step_no = sys.argv[2]

scale_factor = 0

if (step_no == "247"):
    scale_factor = 0.4985

if (step_no == "300"):
    scale_factor = 0.6086

if (step_no == "347"):
    scale_factor = 0.6974

if (step_no == "401"):
    scale_factor = 0.8051

if (step_no == "453"):
    scale_factor = 0.9083

if (step_no == "499"):
    scale_factor = 1.0000

    
bias = np.loadtxt('STEP'+step_no+'/b1_smooth_sod_test.'+step_no+'.txt')

r_min = 0

if int(bin_no) < 7:
    for model_no in range(11,37):
        if step_no == "499":
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L2100_sod_{:5.4f}_{:4.2f}.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        else:
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L2100_sod_{:5.4f}_{:4.2f}_test.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        r_2100,xi_2100,count = np.loadtxt(filename).T
        cutoff = np.where(r_2100 > 0) # want the first non-zero r value
        if (r_min < r_2100[cutoff[0][0]]):
            r_min = r_2100[cutoff[0][0]]





else:
    for model_no in range(11,37):
        if step_no == "499":
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L5000_sod_{:5.4f}_{:4.2f}.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        else:
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L5000_sod_{:5.4f}_{:4.2f}_test.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        r_5000,xi_5000,count = np.loadtxt(filename).T
        cutoff = np.where(r_5000 > 0) # want the first non-zero r value
        if (r_min < r_5000[cutoff[0][0]]):
            r_min = r_5000[cutoff[0][0]]


r_min = np.round(r_min+0.01,2)

print(r_min)




xi = np.zeros([nmodels, nbins])
r = np.logspace(np.log10(r_min), np.log10(100.0), nbins)

#L2100: use bins 0-6 
#L5000: use bins 9-15

if int(bin_no) < 7:
    for model_no in range(30,31):
        fftw_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L2100_{:5.4f}_fftw_test.dat'.format(model_no, model_no, scale_factor)

        r_m_fftw_L2100, xi_m_fftw_L2100 = np.loadtxt(fftw_filename).T

        fftw_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L5000_2048_{:5.4f}_fftw_test.dat'.format(model_no, model_no, scale_factor)

        r_m_fftw_L5000, xi_m_fftw_L5000 = np.loadtxt(fftw_filename).T

        
        corrfunc_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L2100_{:5.4f}_sub_test.dat'.format(model_no, model_no, scale_factor)

        r_m_corrfunc, xi_m_corrfunc, count = np.loadtxt(corrfunc_filename).T
        
        xi_L2100_interp = interp1d(r_m_fftw_L2100, xi_m_fftw_L2100,'cubic')
        xi_L5000_interp = interp1d(r_m_fftw_L5000, xi_m_fftw_L5000,'cubic')

        
        b = bias[model_no-11,int(bin_no)]
        if step_no == "499":
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L2100_sod_{:5.4f}_{:4.2f}.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        else:
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L2100_sod_{:5.4f}_{:4.2f}_test.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        r_2100,xi_2100,count = np.loadtxt(filename).T

       
        f = interp1d(r_2100,r_2100*r_2100*xi_2100/b/b)

        if (int(bin_no) == 0):
            xi_smooth = savgol_filter(f(r), 5, 3)
        elif (int(bin_no) < 4):
            xi_smooth = savgol_filter(f(r), 7, 3)
        else:
            xi_smooth = savgol_filter(f(r), 11, 3)
            
        plt.plot(r_2100[4:], np.sqrt(xi_2100[4:]/xi_L2100_interp(r_2100[4:])), label='fftw (L2100)')
        plt.plot(r_2100[7:], np.sqrt(xi_2100[7:]/xi_L5000_interp(r_2100[7:])), label='fftw (L5000)')
        plt.plot(r_2100, np.sqrt(xi_2100/xi_m_corrfunc), label='corrfunc')
        plt.plot(r_2100, b*np.ones(len(r_2100)), label='linear bias')

        
#        xmax = (np.max(np.where(r_2100 < 100)))

        xi[model_no-11,] = xi_smooth

else:
    for model_no in range(30,31):
        fftw_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L5000_2048_{:5.4f}_fftw_test.dat'.format(model_no, model_no, scale_factor)

        r_m_fftw_L5000, xi_m_fftw_L5000 = np.loadtxt(fftw_filename).T

        fftw_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L2100_{:5.4f}_fftw_test.dat'.format(model_no, model_no, scale_factor)

        r_m_fftw_L2100, xi_m_fftw_L2100 = np.loadtxt(fftw_filename).T
        
        corrfunc_filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/xi_M{:0>3}_L2100_{:5.4f}_sub_test.dat'.format(model_no, model_no, scale_factor)

        r_m_corrfunc, xi_m_corrfunc, count = np.loadtxt(corrfunc_filename).T
        
        xi_L2100_interp = interp1d(r_m_fftw_L2100, xi_m_fftw_L2100,'cubic')
        xi_L5000_interp = interp1d(r_m_fftw_L5000, xi_m_fftw_L5000,'cubic')
        
        b = bias[model_no-11,int(bin_no)]
        if step_no == "499":
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L5000_sod_{:5.4f}_{:4.2f}.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        else:
            filename = '/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{}/xi_M{:0>3}_L5000_sod_{:5.4f}_{:4.2f}_test.dat'.format(model_no, step_no, model_no, scale_factor, mass[int(bin_no)])
        r_5000,xi_5000,count = np.loadtxt(filename).T

        f = interp1d(r_5000,r_5000*r_5000*xi_5000/b/b)

        xi_smooth = savgol_filter(f(r), 11, 3) # used to be 9,3
        if (int(bin_no)==15):
            xi_smooth = savgol_filter(f(r), 15, 3)

#        plt.plot(r,(xi_smooth), label='M{:0>3}'.format(model_no))

        plt.plot(r_5000[7:], np.sqrt(xi_5000[7:]/xi_L2100_interp(r_5000[7:])), label='fftw (L2100)')
        plt.plot(r_5000[7:], np.sqrt(xi_5000[7:]/xi_L5000_interp(r_5000[7:])), label='fftw (L5000)')
        plt.plot(r_5000, np.sqrt(xi_5000/xi_m_corrfunc), label='corrfunc')
        plt.plot(r_5000, b*np.ones(len(r_5000)), label='linear bias')
            

        
        xi[model_no-11,] = xi_smooth


plt.legend(loc='best')
plt.xlim((2,150))
#plt.ylim((1.0,3.0))

#np.savetxt('STEP{}/r.sod.{}.{:4.2f}_test.txt'.format(step_no, step_no, mass[int(bin_no)]),r)
#np.savetxt('STEP{}/xi.sod.{}.{:4.2f}_test.txt'.format(step_no, step_no, mass[int(bin_no)]),xi)

#plt.axis([1.,200.,1e-4,25.])
plt.show()
    
