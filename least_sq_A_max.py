#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130805_Hopkins_y_with_ymax_inset.pdf')
x_nopeak,y_nopeak = numpy.loadtxt('Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_peak,y_peak = numpy.loadtxt('Y24903L_pk.PRN', unpack=True, usecols = [0,1])# L => eps2

z_y = linspace(0,50,500)
sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))
yA_dep_exp = numpy.empty(len(z_y))
norm_yA_dep = numpy.empty(len(z_y))
A_dep = linspace(10,11,12)#[2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
#A_dep = [2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
#A_dep = [11.4*1.4,11.4*1.45,11.4*1.455,11.4*1.5]#,2.,4.,5.]#linspace(5,36,6)
#A_dep = [10.8*1.4,10.8*1.45,10.8*1.455,10.8*1.5]#,2.,4.,5.]#linspace(5,36,6)
#A_dep = [8.4*1.4, 8.4*1.45, 8.4*1.5]#,2.,4.,5.]#linspace(5,36,6)
y_0 = numpy.empty(len(A_dep))
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0
max_yA = numpy.empty(len(A_dep))
max_xA = numpy.empty(len(A_dep))
total = numpy.empty(len(A_dep))

eiz_npk = numpy.loadtxt('130710_eiz_output_eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_pek = numpy.loadtxt('130710_eiz_output_eiz_peak.txt', unpack=True, usecols = [0])
n = arange(0,500)
deiz = eiz_pek - eiz_npk
de_e = deiz/eiz_pek
#pl.figure()
#pl.plot(n, de_e)
#pl.show()
Aiz = numpy.empty(len(eiz_npk))

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
colors = ['b','g','r','c','m','y','k']
for r in range(len(A_dep)):
    labels=r'$A\, =\, %.2f$' % A_dep[r]
    for t in range(len(z_y)):
        sum_Li2_A_dep[t] = 0.0
        for u in arange(1,101):
            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
        yA_dep[t] = (sum_Li2_A_dep[t])*\
        ((1. - ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2))))**2)/\
        ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2)))))*\
        (1./(1.+z_y[t]**2))/(1. + A_dep[r]/(1+z_y[t]**2))

        yA_dep_exp[t] = (sum_Li2_A_dep[t])*\
        ((1. - ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2))))**2)/\
        ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2)))))#*\
        #(1./(1.+z_y[t]**2))/(1. + A_dep[r]/(1+z_y[t]**2))
	Aiz = [2.0*y_dep*dee for y_dep,dee in zip(yA_dep_exp,de_e)]
    ax.plot(n/12.,Aiz,linestyle = ':',color = 'y')#,  color = colors[r], label = labels)
    total = trapz(yA_dep,z_y)
    print total
    ax.plot(z_y,yA_dep)#,  color = colors[r], label = labels)
pl.xlabel(r'$\xi/\omega_{0}$', size = 21)
pl.ylabel(r'$ A(\xi/\omega_{0} ) $', size = 21)


pl.show()
pl.close()        


