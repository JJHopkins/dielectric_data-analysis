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

#x_nopeak = numpy.loadtxt('data/x_nopeak.txt', unpack=True, usecols = [0])
#y_nopeak = numpy.loadtxt('data/a-sio2-eps2.txt', unpack=True, usecols = [1])
#
#x_peak = numpy.loadtxt('data/x_peak.txt', unpack=True, usecols = [0])
#y_peak = numpy.loadtxt('data/Y24204L.txt', unpack=True, usecols = [1])# L => eps2
#
#
#eiz_nopeak = numpy.loadtxt('data/eiz_nopeak.txt', unpack=True, usecols = [0])
#eiz_peak = numpy.loadtxt('data/eiz_peak.txt', unpack=True, usecols = [0])
#
#
#x_ymax, y_ymax = numpy.loadtxt('output/130729_Hopkins_ymax_output.txt', unpack=True, usecols = [0,1])

#coeff = (2.41*1e14) # in rad/s
#n = arange(0,250)
#z = n * coeff
#A = 4.0 #[1,4,6]
#B = 7.0 #[1,4,6]
#C = 10.0 #[1,4,6]
#sum_Li2 = numpy.empty(len(z))
#dF = numpy.empty(len(z))
#sum_Li2_yA = numpy.empty(len(z_y))
#sum_Li2_yB = numpy.empty(len(z_y))
#sum_Li2_yC = numpy.empty(len(z_y))

#yA = numpy.empty(len(z_y))
#yB = numpy.empty(len(z_y))
#yC = numpy.empty(len(z_y))
#A_dep = linspace(1.,11.,300)

#sum_Li3 = numpy.empty(len(z))
#sum_Li3_p = numpy.empty(len(z))
#F_3 = numpy.empty(len(z))
#F_3_p = numpy.empty(len(z))

#z_y = linspace(0,10,250)
z_y = linspace(0,50,500)
sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))
yA_dep_exp = numpy.empty(len(z_y))
norm_yA_dep = numpy.empty(len(z_y))
A_dep = linspace(0,20,50)#[2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
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

eiz_npk = numpy.loadtxt('130710_eiz_output_eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_pek = numpy.loadtxt('130710_eiz_output_eiz_peak.txt', unpack=True, usecols = [0])
n = arange(0,500)
#deiz = eiz_pek/eiz_pek[0] - eiz_npk/eiz_npk[0]
deiz = eiz_pek - eiz_npk
#de_e = deiz/(eiz_npk/eiz_npk[0])
de_e = deiz/eiz_pek
#pl.figure()
#pl.plot(n, de_e)
#pl.show()
Aiz = numpy.empty(len(eiz_npk))
#for p in range(len(z_y)):
#    sum_Li2_yA[p] = 0.0
#    sum_Li2_yB[p] = 0.0
#    sum_Li2_yC[p] = 0.0
#    for q in arange(1,101):
#        #sum_Li2_yA[p] += (((A/(1+z_[p]))**q)/(q**2)# +  (((a/(a+2.*z_y[p]**2 +2)))**q)/(q**2))
#        sum_Li2_yA[p] += (((A/(A+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#        sum_Li2_yB[p] += (((B/(B+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#        sum_Li2_yC[p] += (((C/(C+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#    yA[p] = sum_Li2_yA[p]*((A/(A+2.*z_y[p]**2 +2))**(-1)*(1 - (A/(A+2.*z_y[p]**2 +2))**2))   
#    yB[p] = sum_Li2_yB[p]*((B/(B+2.*z_y[p]**2 +2))**(-1)*(1 - (B/(B+2.*z_y[p]**2 +2))**2))  
#    yC[p] = sum_Li2_yC[p]*((C/(C+2.*z_y[p]**2 +2))**(-1)*(1 - (C/(C+2.*z_y[p]**2 +2))**2)) 


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
	Aiz = [y_dep*dee for y_dep,dee in zip(yA_dep_exp,de_e)]
#        yA_dep[t] = (sum_Li2_A_dep[t]*\
#        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
#        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

#    max_xA[r]= z_y[yA_dep.argmax()]
    
    # MAIN PLOT
    # ----------------------------------------------------------
    #ax.plot(z_y,yA_dep/yA_dep[0],  color = colors[r], label = labels)
    #ax.plot(n/10.4,Aiz/Aiz[0],linestyle = ':',color = 'k')#,  color = colors[r], label = labels)
    #ax.plot(n/ 11.0,Aiz,linestyle = ':',color = 'b')#,  color = colors[r], label = labels)
    #ax.plot(n/(8.4*3.14/2),Aiz,linestyle = ':',color = 'g')#,  color = colors[r], label = labels)
    #ax.plot(n/12.0,Aiz/Aiz[0],linestyle = ':',color = 'r')#,  color = colors[r], label = labels)
    ax.plot(n/8.4,Aiz/Aiz[0],linestyle = ':',color = 'y')#,  color = colors[r], label = labels)

#        yA_dep[t] = (sum_Li2_A_dep[t]*\
#        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
#        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

#    max_xA[r]= z_y[yA_dep.argmax()]
    
    # MAIN PLOT
    # ----------------------------------------------------------
    #ax.plot(z_y,yA_dep/yA_dep[0],  color = colors[r], label = labels)
    #ax.plot(z_y,yA_dep/yA_dep[0])#,  color = colors[r], label = labels)
    ax.plot(z_y,yA_dep/yA_dep[0])#,  color = colors[r], label = labels)
    #ax.plot(z_y,yA_dep)#/yA_dep[0])#,  color = colors[r], label = labels)



pl.xlabel(r'$\xi/\omega_{0}$', size = 21)
pl.ylabel(r'$ A(\xi/\omega_{0} ) $', size = 21)

pl.savefig('many_max_A.png')

#####eiz_npk = numpy.loadtxt('130710_eiz_output_eiz_nopeak.txt', unpack=True, usecols = [0])
#####eiz_pek = numpy.loadtxt('130710_eiz_output_eiz_peak.txt', unpack=True, usecols = [0])
#####n = arange(0,500)
######deiz = eiz_pek/eiz_pek[0] - eiz_npk/eiz_npk[0]
#####deiz = eiz_pek - eiz_npk
######de_e = deiz/(eiz_npk/eiz_npk[0])
#####de_e = deiz/eiz_npk
#####pl.figure()
#####pl.plot(n, de_e)
#####pl.show()
#####Aiz = numpy.empty(len(eiz_npk))
#####
#####pl.figure()
#####for r in range(len(A_dep)):
#####    labels=r'$A\, =\, %.2f$' % A_dep[r]
#####    for t in range(len(z_y)):
#####        sum_Li2_A_dep[t] = 0.0
#####        for u in arange(1,101):
#####            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
#####        yA_dep[t] = (sum_Li2_A_dep[t])*\
#####        ((1. - ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2))))**2)/\
#####        ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2)))))#*\
#####        #(1./(1.+z_y[t]**2))/(1. + A_dep[r]/(1+z_y[t]**2))
#####	Aiz = [y_dep*dee for y_dep,dee in zip(yA_dep,de_e)]
######        yA_dep[t] = (sum_Li2_A_dep[t]*\
######        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
######        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\
#####
######    max_xA[r]= z_y[yA_dep.argmax()]
#####    
#####    # MAIN PLOT
#####    # ----------------------------------------------------------
#####    #ax.plot(z_y,yA_dep/yA_dep[0],  color = colors[r], label = labels)
#####    ax.plot(,Aiz)#,  color = colors[r], label = labels)
#####pl.show()




## INSET
## ----------------------------------------------------------
#ax_inset = fig.add_axes([0.5,0.5,0.36,0.36])
##ax_inset = fig.add_axes([0.55,0.55,0.33,0.33])
## NB, if you change the above eqn for y, you must regenerate this input:
#ax_inset.plot(x_ymax,y_ymax, color = 'k', linestyle = '-')
#pl.xlabel(r'$\frac{2C}{\pi}$',size = 14)#\,=\,\frac{2}{\pi}C
#pl.ylabel(r'$\rm{MAX}$ $[A(\xi/\omega_{0})] $',size = 14)
#pl.tick_params(labelsize = 'small')
##ax.legend(loc = 'best')
#pl.savefig('fontsize-test.jpg')
#
#pp.savefig()
pl.show()
#pp.close()

pl.close()        

