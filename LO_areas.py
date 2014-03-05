#!/usr/bin/python
"""
    A simple example using scipy curve_fit to fit data from a file.

    The example provided is a fit of Gaussian or Lorentzian functions
    
    To fit your own data, you need to change:
    (1) def func(x,*p) to return the function you are trying to fit,
    (2) the name of the data file read in by numpy.loadtxt,
    (3) the initial p0 values in the scipy.optimize.curve_fit call.
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
"""

import matplotlib.pyplot as plt
from matplotlib import axis as ax
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib               
import numpy as np                    # http://numpy.scipy.org/
import scipy                    # http://scipy.org/
import scipy.optimize, scipy.special, scipy.stats
from matplotlib import pyplot as pl
import scipy.integrate as grt
import pylab

from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130814_Hopkins_contour_xiMax_of_r_A.pdf')

def lorentzian(x,height,center,fullwidth) :
    # A lorentzian peak with:
    #   Peak height above background : p[0]
    #   Central value                : p[1]
    #   Full Width at Half Maximum   : p[2]
    #return (p[0]/numpy.pi)/(1.0+((x-p[1])/p[2])**2)
    return (height)/(1.0+((x-center)/fullwidth)**2)

## 1st peak
#xMin0 = 0.0
#xMax0 = 10.5
#p0_guess = (3.91360,10.3022,0.40694)
#
#y_fitData=y_data_expt[(x_data_expt>xMin0)*(x_data_expt<xMax0)]
#x_fitData=x_data_expt[(x_data_expt>xMin0)*(x_data_expt<xMax0)]
#
#exptPeak0Params, exptPeak0CovMat = scipy.optimize.curve_fit(lorentzian,
#        x_fitData, y_fitData, p0=p0_guess)
#
#exptPeak0Lorentz=lorentzian(x_data_expt,*exptPeak0Params)
hs = np.linspace(0.5,5.0,6)
w_0 = 50.0
fwhms = np.linspace(0.05,1.5,6)

xs = np.linspace(1,100,10000)
y = np.empty(len(xs))
area = np.empty(len(hs))

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
colors = ['b','g','r','c','m','y','k']
#lines = [':','-.','--',':','-.','--']
lines = ['--','-','--','-','--','-']
height=np.zeros(shape=(len(hs),len(fwhms)))
for i,h in enumerate(hs):
    for j,f in enumerate(fwhms):
        y = [lorentzian(x,h,w_0,f) for x in xs]
    #area = [np.trapz(y,x) for x in xs]
        area = np.trapz(y,xs)# for x in xs]
        height[i,j] = area#z_y[yA_dep==max(yA_dep)] 
        print h,f,area
        pl.plot(xs,y, color = colors[i], linestyle = lines[j])
    #pl.figure()
pl.axis([40.0,60.0,0.0,5.5])
pl.title(r'$LO peaks$')#$Area under Lorentzian Oscillator peak$')
pl.savefig('plots/1301219_LOs.png', dpi = 300)
#pl.plot(hs,area)
pl.show()


X,Y = np.meshgrid(hs,fwhms)
pl.figure()
#levels = [0.20,0.40,0.60]#,0.80,1.00]
#contourf(X,Y,height, 100, extend = 'both', cmap = hot());cbar = colorbar(); clim(0,1.2);

#CS=contour(X,Y,height, levels)#, colors = 'k')
#clabel(CS,fmt = '%2.1f',colors = 'k',fontsize = 18)
#xlabel(r'$A$', size = 22)#\,=\,\frac{2}{\pi}C$', size = 20)
#ylabel(r'$r$', size = 22)
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cbar.add_lines(CS)

pylab.contourf(X,Y,height,1000)#,cmap = cm.hot())
CS = pylab.contour(X,Y,height, colors = 'k')
man_loc = [(1,1),(2,2),(3,3),(4,4)]
pl.clabel(CS, inline =1,fmt = '%2.1f', fontsize = 18,color = 'k')#, manual = man_loc)
pl.title(r'$Area under Lorentzian Oscillator peak$')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$0<\frac{\xi}{\omega_0}<2$', size = 'x-large')
#pl.ylabel(r'$y(\xi)$', size = 'x-large')# \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
pl.xlabel(r'$height of peak$', size = 21)#\frac{2}{\pi}C$', size = 21)
pl.ylabel(r'$FWHM$', size = 21)
pl.savefig('plots/1301219_LO_areas_contour.png', dpi = 300)
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
##figure()
##contourf(X,Y,h, 1000, cmap = hot())
#surf = ax.plot_surface(X,Y,height)#, rstride = 1, cstride = 1,alpha = 0.9)#, cmap = cm.hot, linewidth = 0.01, antialiased = True, shade = False)# True)#, cmap = hot()
#ax.set_zlim(-1.01, 1.01)
#
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
#fig.colorbar(surf, shrink=0.5, aspect=5)
##ax.plot_surface(X,Y,height)#,10, rstride = 6, cstride = 6,alpha = 1.0, cmap = cm.hot, linewidth = 0.001, antialiased = True, shade = False)# True)#, cmap = hot()
#fig.colorbar(surf)
#savefig('plots/130815_Hopkins_contour_deltaA_r_xiw0.png', dpi = 300)
pl.show()
pl.close()  
       
