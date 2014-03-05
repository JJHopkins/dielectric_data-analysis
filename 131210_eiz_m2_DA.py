#!/usr/bin/python
#import matplotlib               
#from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from scipy.optimize import leastsq
from random import *

pl.close('all')
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
ns = arange(0,500)
z = ns * coeff

x_nopeak  ,y_nopeak  = loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_042_peak,y_042_peak=loadtxt('data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak,y_062_peak=loadtxt('data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak,y_082_peak=loadtxt('data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak,y_102_peak=loadtxt('data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak,y_132_peak=loadtxt('data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak,y_182_peak=loadtxt('data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2

noise0 = normal(0,0.10,len(x_nopeak))  
noise1 = normal(0,0.10,len(x_042_peak))
noise2 = normal(0,0.10,len(x_062_peak))
noise3 = normal(0,0.10,len(x_082_peak))
noise4 = normal(0,0.10,len(x_102_peak))
noise5 = normal(0,0.10,len(x_132_peak))
noise6 = normal(0,0.10,len(x_182_peak))

y_npk_meas = [y_no  + nois0 for y_no ,nois0 in zip(y_nopeak  ,noise0)]# normal(0,0.02,len(x_nopeak))  
y_042_meas = [y_042 + nois1 for y_042,nois1 in zip(y_042_peak,noise1)]# normal(0,0.02,len(x_042_peak))
y_062_meas = [y_062 + nois2 for y_062,nois2 in zip(y_062_peak,noise2)]# normal(0,0.02,len(x_062_peak))
y_082_meas = [y_082 + nois3 for y_082,nois3 in zip(y_082_peak,noise3)]# normal(0,0.02,len(x_082_peak))
y_102_meas = [y_102 + nois4 for y_102,nois4 in zip(y_102_peak,noise4)]# normal(0,0.02,len(x_102_peak))
y_132_meas = [y_132 + nois5 for y_132,nois5 in zip(y_132_peak,noise5)]# normal(0,0.02,len(x_132_peak))
y_182_meas = [y_182 + nois6 for y_182,nois6 in zip(y_182_peak,noise6)]# normal(0,0.02,len(x_182_peak))


pl.figure()
pl.plot(x_nopeak  , y_npk_meas, label =r'$ab\, initio$' )  
pl.plot(x_042_peak, y_042_meas,color= 'c', label =r'$\omega_{0} = 0.64$')
pl.plot(x_062_peak, y_062_meas,color= 'r', label =r'$\omega_{0} = 0.94$')
pl.plot(x_082_peak, y_082_meas,color= 'g', label =r'$\omega_{0} = 1.25$')
pl.plot(x_102_peak, y_102_meas,color= 'm', label =r'$\omega_{0} = 1.55$')
pl.plot(x_132_peak, y_132_meas,color= 'y', label =r'$\omega_{0} = 2.01$')
pl.plot(x_182_peak, y_182_meas,color= 'k', label =r'$\omega_{0} = 2.76$')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
pl.savefig('plots/eps2_m2_sigma10.png', dpi = 300)
pl.show()

# MAKE EMPTY LISTS:
#----------------------------------------------------------
eiz_nopeak   = empty(len(z))
eiz_042_peak = empty(len(z))
eiz_062_peak = empty(len(z))
eiz_082_peak = empty(len(z))
eiz_102_peak = empty(len(z))
eiz_132_peak = empty(len(z))
eiz_182_peak = empty(len(z))
eizs = empty(len(z))
noise=empty(len(x_nopeak))
meas_ys=empty(len(x_nopeak))

listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_042_peak_arg=empty(len(x_042_peak))
    eiz_062_peak_arg=empty(len(x_062_peak))
    eiz_082_peak_arg=empty(len(x_082_peak))
    eiz_102_peak_arg=empty(len(x_102_peak))
    eiz_132_peak_arg=empty(len(x_132_peak))
    eiz_182_peak_arg=empty(len(x_182_peak))
    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_npk_meas[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for k in range(len(x_042_peak)):
        eiz_042_peak_arg[k]=x_042_peak[k]*y_042_meas[k]/(x_042_peak[k]**2 + z[j]**2)
    eiz_042_peak[j] = 1. + (2./pi)*trapz(eiz_042_peak_arg,x_042_peak)    

    for m in range(len(x_062_peak)):
        eiz_062_peak_arg[m]=x_062_peak[m]*y_062_meas[m]/(x_062_peak[m]**2 + z[j]**2)
    eiz_062_peak[j] = 1. + (2./pi)*trapz(eiz_062_peak_arg,x_062_peak)    

    for n in range(len(x_082_peak)):
        eiz_082_peak_arg[n]=x_082_peak[n]*y_082_meas[n]/(x_082_peak[n]**2 + z[j]**2)
    eiz_082_peak[j] = 1. + (2./pi)*trapz(eiz_082_peak_arg,x_082_peak)    

    for p in range(len(x_102_peak)):
        eiz_102_peak_arg[p]=x_102_peak[p]*y_102_meas[p]/(x_102_peak[p]**2 + z[j]**2)
    eiz_102_peak[j] = 1. + (2./pi)*trapz(eiz_102_peak_arg,x_102_peak)    

    for q in range(len(x_132_peak)):
        eiz_132_peak_arg[q]=x_132_peak[q]*y_132_meas[q]/(x_132_peak[q]**2 + z[j]**2)
    eiz_132_peak[j] = 1. + (2./pi)*trapz(eiz_132_peak_arg,x_132_peak)    

    for s in range(len(x_182_peak)):
        eiz_182_peak_arg[s]=x_182_peak[s]*y_182_meas[s]/(x_182_peak[s]**2 + z[j]**2)
    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    
#         
savetxt("data/eiz_npk_meas2.txt", eiz_nopeak)
savetxt("data/eiz_042_meas2.txt", eiz_042_peak)
savetxt("data/eiz_062_meas2.txt", eiz_062_peak)
savetxt("data/eiz_082_meas2.txt", eiz_082_peak)
savetxt("data/eiz_102_meas2.txt", eiz_102_peak)
savetxt("data/eiz_132_meas2.txt", eiz_132_peak)
savetxt("data/eiz_182_meas2.txt", eiz_182_peak)
#
pl.close()





