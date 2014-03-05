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
pl.savefig('plots/eps2_m0_sigma10.png', dpi = 300)
pl.show()







##ys=zeros(shape=(len(y_nopeak),len(7)))
#ys = [y_nopeak   + normal(0,0.05,len(x_nopeak))  
#,y_042_peak + normal(0,0.05,len(x_042_peak))
#,y_062_peak + normal(0,0.05,len(x_062_peak))
#,y_082_peak + normal(0,0.05,len(x_082_peak))
#,y_102_peak + normal(0,0.05,len(x_102_peak))
#,y_132_peak + normal(0,0.05,len(x_132_peak))
#,y_182_peak + normal(0,0.05,len(x_182_peak))]
#
#xs = [ x_nopeak  
#,x_042_peak
#,x_062_peak
#,x_082_peak
#,x_102_peak
#,x_132_peak
#,x_182_peak]

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
#noise = zeros(shape=(shape(ys))
#ys=zeros(shape=(len(y_nopeak),len(7)))
#pl.figure()
#for a in range(0,7):
#    noise = normal(0,0.05,len(x_nopeak))
#    print noise[100]
#    meas_ys = ys[a]+noise# for noise in noises]
    #pl.plot(x_042_peak,noise)
#    pl.plot(x_042_peak,meas_ys)
#print shape(meas_ys)
#for d,meas_y in enumerate(meas_ys):
#for b in range(len(z)):
#noise = normal(0,0.05,len(xs[a]))
#for a in range(0,7):
#    noise = normal(0,0.05,len(xs[a]))
    #print noise[100]
#    meas_y = ys[a]+noise[a]# for noise in noises]
#    meas_ys = [meas_y[a]]	
#    eizs = empty(len(z))
#    for b in range(len(z)):
#        eiz_arg=empty(len(x_nopeak))
#       # print a,c,b
#        for c in range(len(x_nopeak)):
#        #print meas_ys[1000]
#            eiz_arg[c] = x_nopeak[c]*meas_ys[c] / (x_nopeak[c]**2 + z[b]**2)
#            #print a,b,c
#        eizs[b] = 1. + (2./pi) * trapz(eiz_arg,x_nopeak)
#
#pl.plot(ns,eizs)
#pl.show()

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
#pl.figure()
#pl.plot(ns, eiz_nopeak, label =r'$ab\, initio$')  
#pl.plot(ns, eiz_042_peak, label =r'$\omega_{0} = 0.64$') 
#pl.plot(ns, eiz_062_peak, label =r'$\omega_{0} = 0.94$') 
#pl.plot(ns, eiz_082_peak, label =r'$\omega_{0} = 1.25$') 
#pl.plot(ns, eiz_102_peak, label =r'$\omega_{0} = 1.55$')
#pl.plot(ns, eiz_132_peak, label =r'$\omega_{0} = 2.01$')
#pl.plot(ns, eiz_182_peak, label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$n$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi)\, \rm{for\, ab\, initio\, and\, L.O.\, peaks \,at}\, \omega_{0}$')
#pl.legend(loc = 'best')
#pl.show()
#pl.close()


#for j in range(len(z)):
#    eiz_nopeak_arg=empty(len(x_nopeak))
#    eiz_042_peak_arg=empty(len(x_042_peak))
#    eiz_062_peak_arg=empty(len(x_062_peak))
#    eiz_082_peak_arg=empty(len(x_082_peak))
#    eiz_102_peak_arg=empty(len(x_102_peak))
#    eiz_132_peak_arg=empty(len(x_132_peak))
#    eiz_182_peak_arg=empty(len(x_182_peak))
#    for i in range(len(x_nopeak)):
#        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
#    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)
#
#    for k in range(len(x_042_peak)):
#        eiz_042_peak_arg[k]=x_042_peak[k]*y_042_peak[k]/(x_042_peak[k]**2 + z[j]**2)
#    eiz_042_peak[j] = 1. + (2./pi)*trapz(eiz_042_peak_arg,x_042_peak)    
#
#    for m in range(len(x_062_peak)):
#        eiz_062_peak_arg[m]=x_062_peak[m]*y_062_peak[m]/(x_062_peak[m]**2 + z[j]**2)
#    eiz_062_peak[j] = 1. + (2./pi)*trapz(eiz_062_peak_arg,x_062_peak)    
#
#    for n in range(len(x_082_peak)):
#        eiz_082_peak_arg[n]=x_082_peak[n]*y_082_peak[n]/(x_082_peak[n]**2 + z[j]**2)
#    eiz_082_peak[j] = 1. + (2./pi)*trapz(eiz_082_peak_arg,x_082_peak)    
#
#    for p in range(len(x_102_peak)):
#        eiz_102_peak_arg[p]=x_102_peak[p]*y_102_peak[p]/(x_102_peak[p]**2 + z[j]**2)
#    eiz_102_peak[j] = 1. + (2./pi)*trapz(eiz_102_peak_arg,x_102_peak)    
#
#    for q in range(len(x_132_peak)):
#        eiz_132_peak_arg[q]=x_132_peak[q]*y_132_peak[q]/(x_132_peak[q]**2 + z[j]**2)
#    eiz_132_peak[j] = 1. + (2./pi)*trapz(eiz_132_peak_arg,x_132_peak)    
#
#    for s in range(len(x_182_peak)):
#        eiz_182_peak_arg[s]=x_182_peak[s]*y_182_peak[s]/(x_182_peak[s]**2 + z[j]**2)
#    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    
#
#savetxt("data/x_nopeak.txt", x_nopeak    )
#savetxt("data/x_042_peak.txt", x_042_peak)
#savetxt("data/x_062_peak.txt", x_062_peak)
#savetxt("data/x_082_peak.txt", x_082_peak)
#savetxt("data/x_102_peak.txt", x_102_peak)
#savetxt("data/x_132_peak.txt", x_132_peak)
#savetxt("data/x_182_peak.txt", x_182_peak)
#         
savetxt("data/eiz_npk_meas.txt", eiz_nopeak)
savetxt("data/eiz_042_meas.txt", eiz_042_peak)
savetxt("data/eiz_062_meas.txt", eiz_062_peak)
savetxt("data/eiz_082_meas.txt", eiz_082_peak)
savetxt("data/eiz_102_meas.txt", eiz_102_peak)
savetxt("data/eiz_132_meas.txt", eiz_132_peak)
savetxt("data/eiz_182_meas.txt", eiz_182_peak)
#
#pl.figure()
#pl.plot(  x_nopeak, y_nopeak  , label =r'$ab\, initio$' )  
#pl.plot(x_042_peak, y_042_peak,color= 'c', label =r'$\omega_{0} = 0.64$')
#pl.plot(x_062_peak, y_062_peak,color= 'r', label =r'$\omega_{0} = 0.94$')
#pl.plot(x_082_peak, y_082_peak,color= 'g', label =r'$\omega_{0} = 1.25$')
#pl.plot(x_102_peak, y_102_peak,color= 'm', label =r'$\omega_{0} = 1.55$')
#pl.plot(x_132_peak, y_132_peak,color= 'y', label =r'$\omega_{0} = 2.01$')
#pl.plot(x_182_peak, y_182_peak,color= 'k', label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
#pl.savefig('plots/131210_eps2_eV.png', dpi = 300)
#
#pl.show()
pl.close()



