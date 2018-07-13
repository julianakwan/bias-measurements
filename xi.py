import numpy as np
import treecorr
import sys
from scipy.signal import savgol_filter
import Corrfunc
from Corrfunc.io import read_catalog
from Corrfunc.theory.xi import xi

def read_halo_catalogue(model_no, step_no, bin_no):
    mass_cut = np.arange(13.0,15.0,0.1)
    filename = "/Users/jkwan/emulators/TitanU/xi/M{:0>3d}/STEP{:3d}/halo_catalogue_M{:0>3d}_L2100_sod_0.4985_{:4.2f}.dat".format(int(model_no), int(step_no), int(model_no), mass_cut[int(bin_no)])
    x,y,z = np.loadtxt(filename).T

    print(filename)
#    x,y,z,npart = np.loadtxt('/Users/jkwan/emulators/TitanU/halo_power/M026/halo_catalogue_M026_L2100_fof_1.0000_9.dat').T

    #Only take a subsample
    if (x.size > 1e8):
        index = np.random.randint(x.size, size=int(0.1*x.size))
        x_sub = x[index]
        y_sub = y[index]
        z_sub = z[index]
        x = x_sub
        y = y_sub
        z = z_sub

    
    return x,y,z


def read_matter_particles(model_no, step_no, sub):
    filename="/Users/jkwan/emulators/TitanU/density/M{:0>3}/STEP{:3}/particle_subsample_M{:0>3}_L2100_1.0000_0.dat".format(model_no, step_no, model_no)
    p = np.fromfile(filename,dtype='f')

    p_new = np.reshape(p,(p.size/6,6))
    # there are 6 fields for struct particle_data

    x_m = p_new[:,0] 
    y_m = p_new[:,1]
    z_m = p_new[:,2]

    # the rest are velocities

    
    #Only take a subsample
    if (sub):
#        x_sub = x_m[(x_m < 0.5) & (y_m < 0.5) & (z_m < 0.5)]
#        y_sub = y_m[(x_m < 0.5) & (y_m < 0.5) & (z_m < 0.5)]
#        z_sub = z_m[(x_m < 0.5) & (y_m < 0.5) & (z_m < 0.5)]
        index = np.random.randint(x_m.size, size=int(0.1*x_m.size))
        x_sub = x_m[index]
        y_sub = y_m[index]
        z_sub = z_m[index]
        x_m = x_sub
        y_m = y_sub
        z_m = z_sub

    return x_m, y_m, z_m

def setup_treecorr():
    max_size = 100000

    x_r = np.random.rand(10*max_size) # randoms
    y_r = np.random.rand(10*max_size)
    z_r = np.random.rand(10*max_size)

    cat_r = treecorr.Catalog(x=x_r,y=y_r,z=z_r)
    dd = treecorr.NNCorrelation(nbins=20,min_sep=0.001,max_sep=0.1,bin_slop=0.05)
    dr = treecorr.NNCorrelation(nbins=20,min_sep=0.001,max_sep=0.1,bin_slop=0.05)
    rr = treecorr.NNCorrelation(nbins=20,min_sep=0.001,max_sep=0.1,bin_slop=0.05)

    rr.process(cat_r)
    
    return cat_r, dd, dr, rr

if __name__ == "__main__":

    model_no = sys.argv[1]
    step_no = sys.argv[2]
    nthreads = 2
    jk_no = 2
    
    mass_cut = np.arange(13.0,15.0,0.1)
    rmin = 1.
    rmax = 200.
    rbins = 20
    BoxSize = 2100.

    log_rbins = np.logspace(np.log10(rmin),np.log10(rmax),rbins+1)

    #These are things that stay the same with every run
    cat_r, dd, dr, rr = setup_treecorr()

    nbins = 16
    
    for bin_no in range(4,5):
        x,y,z = read_halo_catalogue(model_no, step_no, bin_no)
#        x,y,z = read_matter_particles(model_no, step_no,False)

#Divide sim into 1/8 cubes
        if jk_no == 0:
            indices = np.where((x < 0.5) & (y < 0.5) & (z < 0.5))
            x_sub = x[indices[0]]
            y_sub = y[indices[0]]
            z_sub = z[indices[0]]

        if jk_no == 1:
            indices = np.where((x < 0.5) & (y < 0.5) & (z > 0.5))
            x_sub = x[indices[0]]
            y_sub = y[indices[0]]
            z_sub = z[indices[0]]-0.5         # Recenter the subcube. This has to change depending on the subsample taken

        if jk_no == 2:
            indices = np.where((x < 0.5) & (y > 0.5) & (z < 0.5))
            x_sub = x[indices[0]]
            y_sub = y[indices[0]]-0.5        
            z_sub = z[indices[0]]

        if jk_no == 3:
            indices = np.where((x < 0.5) & (y > 0.5) & (z > 0.5))
            x_sub = x[indices[0]]
            y_sub = y[indices[0]]-0.5
            z_sub = z[indices[0]]-0.5

        if jk_no == 4:
            indices = np.where((x > 0.5) & (y < 0.5) & (z < 0.5))
            x_sub = x[indices[0]]-0.5
            y_sub = y[indices[0]]
            z_sub = z[indices[0]]

        if jk_no == 5:
            indices = np.where((x > 0.5) & (y < 0.5) & (z > 0.5))
            x_sub = x[indices[0]]-0.5
            y_sub = y[indices[0]]
            z_sub = z[indices[0]]-0.5

        if jk_no == 6:
            indices = np.where((x > 0.5) & (y > 0.5) & (z < 0.5))
            x_sub = x[indices[0]]-0.5
            y_sub = y[indices[0]]-0.5
            z_sub = z[indices[0]]

        if jk_no == 7:
            indices = np.where((x > 0.5) & (y > 0.5) & (z > 0.5))
            x_sub = x[indices[0]]-0.5
            y_sub = y[indices[0]]-0.5
            z_sub = z[indices[0]]-0.5
            



        x_sub *= BoxSize
        y_sub *= BoxSize
        z_sub *= BoxSize


        
#        cat = treecorr.Catalog(x=x,y=y,z=z)
        #units are arcmin by default, but pretend MPc

#        dd.process(cat)
#        dr.process(cat,cat_r)


#        xi,varxi = dd.calculateXi(rr,dr)

 #       xi_smooth = savgol_filter(xi, 11, 5)
        
 #       r = np.exp(dd.logr)


 
        results = xi(0.5*BoxSize, nthreads, log_rbins, x_sub,y_sub,z_sub,output_ravg=True)

        output_filename = "/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{:3}/xi_M{:0>3}_L2100_sod_0.4985_{:4.2f}_test.{}.dat".format(model_no, step_no, model_no, mass_cut[int(bin_no)],jk_no)
#        output_filename = "/Users/jkwan/emulators/TitanU/xi/M{:0>3}/STEP{:3}/xi_mm_M{:0>3}_L5000_1.0000_test.dat".format(model_no, step_no, model_no)

#        np.savetxt(output_filename, np.column_stack([r*BoxSize, xi]))
        np.savetxt(output_filename, np.column_stack([results['ravg'], results['xi'], results['npairs']]))
                            
