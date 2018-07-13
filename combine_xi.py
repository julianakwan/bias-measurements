import numpy as np
import matplotlib.pyplot as plt

mass_bin = np.concatenate([np.linspace(13.,13.6,7), np.linspace(13.9,14.5,7)])
step_no = [247, 300, 347, 401, 453, 499]
xi = np.zeros([26,6*20])
r  = np.zeros([6*20])

#for i in range(0,len(mass_bin)):
#    xi[:,i*20:(i+1)*20]= np.loadtxt("STEP499/xi.sod.499.{:4.2f}_test.txt".format(mass_bin[i]))
#    r[i*20:(i+1)*20] = np.loadtxt("STEP499/r.sod.499.{:4.2f}_test.txt".format(mass_bin[i]))
#    plt.plot(100*i+r[i*20:(i+1)*20],xi[2,i*20:(i+1)*20])

for i in range(0,1):
    for j in range(0,6):
        xi[:,j*20:(j+1)*20]= np.loadtxt("STEP{:3d}/xi.sod.{:3d}.{:4.2f}_test.txt".format(step_no[j], step_no[j], mass_bin[i]))
        r[j*20:(j+1)*20] = np.loadtxt("STEP{:3d}/r.sod.{:3d}.{:4.2f}_test.txt".format(step_no[j], step_no[j], mass_bin[i]))
        for k in range(0,26):
            plt.plot(100*j+r[j*20:(j+1)*20],xi[k,j*20:(j+1)*20])
    

plt.show()
#np.savetxt("STEP499/xi.sod.499.all.txt",xi)
#np.savetxt("STEP499/r.sod.499.all.txt",r)    

np.savetxt("xi.sod.{:4.2f}.all.txt".format(mass_bin[i]), xi)
np.savetxt("r.sod.{:4.2f}.all.dat".format(mass_bin[i]), r)
