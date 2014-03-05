#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from scipy.optimize import leastsq
from random import *
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130729_Hopkins_shifted peaks_eiz.pdf')

pl.close('all')
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,100)
z = n * coeff

x_nopeak_eV, y_nopeak = loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_042_peak_eV,y_042_peak=loadtxt('data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=loadtxt('data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=loadtxt('data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=loadtxt('data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=loadtxt('data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=loadtxt('data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2

#ys=zeros(shape=(len(y_nopeak),len(7)))

# MAKE EMPTY LISTS:
#----------------------------------------------------------
x_nopeak  = empty(len(x_nopeak_eV))
x_042_peak= empty(len(x_042_peak_eV))
x_062_peak= empty(len(x_062_peak_eV))
x_082_peak= empty(len(x_082_peak_eV))
x_102_peak= empty(len(x_102_peak_eV))
x_132_peak= empty(len(x_132_peak_eV))
x_182_peak= empty(len(x_182_peak_eV))

eiz_nopeak   = empty(len(z))
eiz_042_peak = empty(len(z))
eiz_062_peak = empty(len(z))
eiz_082_peak = empty(len(z))
eiz_102_peak = empty(len(z))
eiz_132_peak = empty(len(z))
eiz_182_peak = empty(len(z))

diff_eiz_042 =  empty(len(z))
diff_eiz_062 =  empty(len(z))
diff_eiz_082 =  empty(len(z))
diff_eiz_102 =  empty(len(z))
diff_eiz_132 =  empty(len(z))
diff_eiz_182 =  empty(len(z))

# Add some noise
noise = normal(0,0.05,len(x_nopeak_eV))
measured_signal = y_042_peak+noise
pl.figure()
pl.plot(x_042_peak_eV,y_042_peak)
pl.plot(x_042_peak_eV,measured_signal)
pl.show()

# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
x_nopeak  = x_nopeak_eV  #*1.5187*10**(15)
x_042_peak= x_042_peak_eV#*1.5187*10**(15)
x_062_peak= x_062_peak_eV#*1.5187*10**(15)
x_082_peak= x_082_peak_eV#*1.5187*10**(15)
x_102_peak= x_102_peak_eV#*1.5187*10**(15)
x_132_peak= x_132_peak_eV#*1.5187*10**(15)
x_182_peak= x_182_peak_eV#*1.5187*10**(15)


listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

#pl.figure()
for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_042_peak_arg=empty(len(x_042_peak))
    eiz_062_peak_arg=empty(len(x_062_peak))
    eiz_082_peak_arg=empty(len(x_082_peak))
    eiz_102_peak_arg=empty(len(x_102_peak))
    eiz_132_peak_arg=empty(len(x_132_peak))
    eiz_182_peak_arg=empty(len(x_182_peak))
    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for k in range(len(x_042_peak)):
        eiz_042_peak_arg[k]=x_042_peak[k]*y_042_peak[k]/(x_042_peak[k]**2 + z[j]**2)
        #print I_w_ew_042
        #pl.figure()
        #pl.plot(x_042_peak, I_w_ew_042)
        #pl.show()
        #pl.close()
    eiz_042_peak[j] = 1. + (2./pi)*trapz(eiz_042_peak_arg,x_042_peak)    

    for m in range(len(x_062_peak)):
        eiz_062_peak_arg[m]=x_062_peak[m]*y_062_peak[m]/(x_062_peak[m]**2 + z[j]**2)
    eiz_062_peak[j] = 1. + (2./pi)*trapz(eiz_062_peak_arg,x_062_peak)    

    for n in range(len(x_082_peak)):
        eiz_082_peak_arg[n]=x_082_peak[n]*y_082_peak[n]/(x_082_peak[n]**2 + z[j]**2)
    eiz_082_peak[j] = 1. + (2./pi)*trapz(eiz_082_peak_arg,x_082_peak)    

    for p in range(len(x_102_peak)):
        eiz_102_peak_arg[p]=x_102_peak[p]*y_102_peak[p]/(x_102_peak[p]**2 + z[j]**2)
    eiz_102_peak[j] = 1. + (2./pi)*trapz(eiz_102_peak_arg,x_102_peak)    

    for q in range(len(x_132_peak)):
        eiz_132_peak_arg[q]=x_132_peak[q]*y_132_peak[q]/(x_132_peak[q]**2 + z[j]**2)
    eiz_132_peak[j] = 1. + (2./pi)*trapz(eiz_132_peak_arg,x_132_peak)    

    for s in range(len(x_182_peak)):
        eiz_182_peak_arg[s]=x_182_peak[s]*y_182_peak[s]/(x_182_peak[s]**2 + z[j]**2)
    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    


    I_w_ew_npk = trapz(x_nopeak*1.5187*10**(-1)*y_nopeak    ) 
    I_w_ew_042 = trapz(x_042_peak*1.5187*10**(-1)*y_042_peak)
    I_w_ew_062 = trapz(x_062_peak*1.5187*10**(-1)*y_062_peak)
    I_w_ew_082 = trapz(x_082_peak*1.5187*10**(-1)*y_082_peak)
    I_w_ew_102 = trapz(x_102_peak*1.5187*10**(-1)*y_102_peak)
    I_w_ew_132 = trapz(x_132_peak*1.5187*10**(-1)*y_132_peak)
    I_w_ew_182 = trapz(x_182_peak*1.5187*10**(-1)*y_182_peak)
    

#pl.plot(x_nopeak,eiz_nopeak_arg,  label = r'$\rm{ab\,initio}$')
#pl.plot(x_182_peak,eiz_182_peak_arg,  label = r'$\omega_{0}=18.2 eV$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi) \rm{for\, ab\, initio\, and\,} \omega_{0}=18.2 eV\, \rm{data}$')
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
#pl.close()

savetxt("data/x_nopeak.txt", x_nopeak    )
savetxt("data/x_042_peak.txt", x_042_peak)
savetxt("data/x_062_peak.txt", x_062_peak)
savetxt("data/x_082_peak.txt", x_082_peak)
savetxt("data/x_102_peak.txt", x_102_peak)
savetxt("data/x_132_peak.txt", x_132_peak)
savetxt("data/x_182_peak.txt", x_182_peak)
         
savetxt("data/eiz_nopeak.txt", eiz_nopeak)
savetxt("data/eiz_042_peak.txt", eiz_042_peak)
savetxt("data/eiz_062_peak.txt", eiz_062_peak)
savetxt("data/eiz_082_peak.txt", eiz_082_peak)
savetxt("data/eiz_102_peak.txt", eiz_102_peak)
savetxt("data/eiz_132_peak.txt", eiz_132_peak)
savetxt("data/eiz_182_peak.txt", eiz_182_peak)


pl.figure()
pl.plot(  x_nopeak, y_nopeak  , label =r'$ab\, initio$' )  
pl.plot(x_042_peak, y_042_peak,color= 'c', label =r'$\omega_{0} = 0.64$')
pl.plot(x_062_peak, y_062_peak,color= 'r', label =r'$\omega_{0} = 0.94$')
pl.plot(x_082_peak, y_082_peak,color= 'g', label =r'$\omega_{0} = 1.25$')
pl.plot(x_102_peak, y_102_peak,color= 'm', label =r'$\omega_{0} = 1.55$')
pl.plot(x_132_peak, y_132_peak,color= 'y', label =r'$\omega_{0} = 2.01$')
pl.plot(x_182_peak, y_182_peak,color= 'k', label =r'$\omega_{0} = 2.76$')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peaks at $\omega_{0}$')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.legend(loc = 'best')
#pl.axis([0,7,0,45])
#pp.savefig()
#pl.savefig('plots/130814_eV_shifted_peaks_eps2.pdf')
pl.savefig('plots/131210_eps2_eV.png', dpi = 300)

pl.figure()
pl.plot(  x_nopeak*1.5187*10**(15)*10**(-16), y_nopeak  , label =r'$ab\, initio$' )  
pl.plot(x_042_peak*1.5187*10**(15)*10**(-16), y_042_peak,color= 'c', label =r'$\omega_{0} = 0.64$')
pl.plot(x_062_peak*1.5187*10**(15)*10**(-16), y_062_peak,color= 'r', label =r'$\omega_{0} = 0.94$')
pl.plot(x_082_peak*1.5187*10**(15)*10**(-16), y_082_peak,color= 'g', label =r'$\omega_{0} = 1.25$')
pl.plot(x_102_peak*1.5187*10**(15)*10**(-16), y_102_peak,color= 'm', label =r'$\omega_{0} = 1.55$')
pl.plot(x_132_peak*1.5187*10**(15)*10**(-16), y_132_peak,color= 'y', label =r'$\omega_{0} = 2.01$')
pl.plot(x_182_peak*1.5187*10**(15)*10**(-16), y_182_peak,color= 'k', label =r'$\omega_{0} = 2.76$')
pl.xlabel(r'$\omega\,\,\,\,[sec^{-1}]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peaks at $\omega_{0}$')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.legend(loc = 'best')
pl.axis([0,7,0,4.1])
#pp.savefig()
pl.savefig('plots/131210_eps2_inverse_sec.png', dpi = 300)
pl.show()
pl.close()

#pl.figure()
#pl.plot(  x_nopeak, y_nopeak  , label =r'ab initio' )  
#pl.plot(x_182_peak, y_182_peak, label =r'$\omega_{0} = 18.2 eV$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peak at $\omega_{0}$ = 18.2 eV')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()
#
#pl.figure()
#pl.plot(z, eiz_nopeak, label =r'$ab\, initio$')  
#pl.plot(z, eiz_042_peak, label =r'$\omega_{0} = 0.64$') 
#pl.plot(z, eiz_062_peak, label =r'$\omega_{0} = 0.94$') 
#pl.plot(z, eiz_082_peak, label =r'$\omega_{0} = 1.25$') 
#pl.plot(z, eiz_102_peak, label =r'$\omega_{0} = 1.55$')
#pl.plot(z, eiz_132_peak, label =r'$\omega_{0} = 2.01$')
#pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi)\, \rm{for\, ab\, initio\, and\, L.O.\, peaks \,at}\, \omega_{0}$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()
#
#pl.figure()
#pl.plot(z, eiz_nopeak, label =r'ab initio')  
#pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 18.2 eV')#2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi)$ for ab initio and L.O. peak at $\omega_{0}$ = 18.2 eV$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()


### calc difference eps2 and eiz for with and without peak
##-----------------------------------------------------------
#diff_eiz_042 = eiz_042_peak - eiz_nopeak
#diff_eiz_062 = eiz_062_peak - eiz_nopeak
#diff_eiz_082 = eiz_082_peak - eiz_nopeak
#diff_eiz_102 = eiz_102_peak - eiz_nopeak
#diff_eiz_132 = eiz_132_peak - eiz_nopeak
#diff_eiz_182 = eiz_182_peak - eiz_nopeak
#
#pl.figure()
#pl.plot(  z/4.2 , diff_eiz_042, color='g', label =r'$\omega_{0} = 0.64$')
#pl.plot(  z/6.2 , diff_eiz_062, color='r', label =r'$\omega_{0} = 0.94$')
#pl.plot(  z/8.2 , diff_eiz_082, color='c', label =r'$\omega_{0} = 1.25$')
#pl.plot(  z/10.2, diff_eiz_102, color='m', label =r'$\omega_{0} = 1.55$')
#pl.plot(  z/13.2, diff_eiz_132, color='y', label =r'$\omega_{0} = 2.01$')
#pl.plot(  z/18.2, diff_eiz_182, color='k', label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$', size = 'x-large')
#pl.ylabel(r'$ \delta\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'Difference between $\epsilon^{synth}(i\xi)$ and $\epsilon^{ab\,int}(i\xi)$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()
#
#pl.figure()
#pl.plot(    x_s*1.5187*10**(15),y_s,  label = r'$ab\,initio$')
#pl.plot(x_042_s*1.5187*10**(15),y_042_s,  label = r'$\omega_0= 0.64$')
#pl.plot(x_062_s*1.5187*10**(15),y_062_s,  label = r'$\omega_0= 0.94$')
#pl.plot(x_082_s*1.5187*10**(15),y_082_s,  label = r'$\omega_0= 1.25$')
#pl.plot(x_102_s*1.5187*10**(15),y_102_s,  label = r'$\omega_0= 1.55$')
#pl.plot(x_132_s*1.5187*10**(15),y_132_s,  label = r'$\omega_0= 2.01$')
#pl.plot(x_182_s*1.5187*10**(15),y_182_s,  label = r'$\omega_0= 2.76$')
#
#pl.plot(  x_nopeak*1.5187*10**(15),  x_nopeak*1.5187*10**(-1)*y_nopeak  , linestyle = ':')#, label =r'$ab\, initio$'      )  
#pl.plot(x_042_peak*1.5187*10**(15),x_042_peak*1.5187*10**(-1)*y_042_peak, linestyle = ':')#, label =r'$\omega_{0} = 0.64$')
#pl.plot(x_062_peak*1.5187*10**(15),x_062_peak*1.5187*10**(-1)*y_062_peak, linestyle = ':')#, label =r'$\omega_{0} = 0.94$')
#pl.plot(x_082_peak*1.5187*10**(15),x_082_peak*1.5187*10**(-1)*y_082_peak, linestyle = ':')#, label =r'$\omega_{0} = 1.25$')
#pl.plot(x_102_peak*1.5187*10**(15),x_102_peak*1.5187*10**(-1)*y_102_peak, linestyle = ':')#, label =r'$\omega_{0} = 1.55$')
#pl.plot(x_132_peak*1.5187*10**(15),x_132_peak*1.5187*10**(-1)*y_132_peak, linestyle = ':')#, label =r'$\omega_{0} = 2.01$')
#pl.plot(x_182_peak*1.5187*10**(15),x_182_peak*1.5187*10**(-1)*y_182_peak, linestyle = ':')#, label =r'$\omega_{0} = 2.76$')
#
#pl.yscale('log')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$')#, size = 'x-large')
#pl.ylabel(r'Strength per unit volume (45.3 A$^{3}$)')#, size = 'x-large')
#pl.title('Sum Rule')
##pl.axis([0.0,22.0,0.0,6.0])
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
#pl.close()


#pl.figure()
#pl.plot(    x_s,y_s,  label = r'$ab\,initio$')
#pl.plot(x_042_s/4.2 ,y_042_s,  label = r'$\omega_0 = 0.64$')
#pl.plot(x_062_s/6.2 ,y_062_s,  label = r'$\omega_0 = 0.94$')
#pl.plot(x_082_s/8.2 ,y_082_s,  label = r'$\omega_0 = 1.25$')
#pl.plot(x_102_s/10.2,y_102_s,  label = r'$\omega_0 = 1.55$')
#pl.plot(x_132_s/13.2,y_132_s,  label = r'$\omega_0 = 2.01$')
#pl.plot(x_182_s/18.2,y_182_s,  label = r'$\omega_0 = 2.76$')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$')#, size = 'x-large')
#pl.ylabel(r'Strength per unit volume (45.3 A$^{3}$)')#, size = 'x-large')
#pl.title('Sum Rule')
##pl.axis([0.0,22.0,0.0,6.0])
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
#pp.close()









