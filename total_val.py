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
#A_dep = linspace(10,11,12)#[2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
A_dep = linspace(0,20,100)#[2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
#A_dep = [2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
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
	Aiz = [2.0*y_dep*dee for y_dep,dee in zip(yA_dep_exp,de_e)]
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
#    ax.plot(n/12.0,Aiz,linestyle = ':',color = 'r')#,  color = colors[r], label = labels)
    #ax.plot(n/12.,Aiz,linestyle = ':',color = 'y')#,  color = colors[r], label = labels)
    total = trapz(yA_dep/yA_dep[0],z_y)
    #print total
#        yA_dep[t] = (sum_Li2_A_dep[t]*\
#        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
#        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

#    max_xA[r]= z_y[yA_dep.argmax()]
    #pl.plot(A_dep,total)
    # MAIN PLOT
    # ----------------------------------------------------------
    #ax.plot(z_y,yA_dep/yA_dep[0],  color = colors[r], label = labels)
    #ax.plot(z_y,yA_dep/yA_dep[0])#,  color = colors[r], label = labels)
    ax.plot(z_y,yA_dep/yA_dep[0])#,  color = colors[r], label = labels)
    #ax.plot(z_y,yA_dep)#/yA_dep[0])#,  color = colors[r], label = labels)

pl.xlabel(r'$\xi/\omega_{0}$', size = 21)
pl.ylabel(r'$ A(\xi/\omega_{0} ) $', size = 21)

totals = [0.0
,0.844235588588
,0.902052850418
,0.959171922439
,1.01581508011
,1.07214213658
,1.12827223495
,1.18429669819
,1.24028706962
,1.2963003875
,1.35238278039
,1.40857199673
,1.46489923268
,1.52139048298
,1.57806755879
,1.63494886788
,1.69205002126
,1.74938431132
,1.80696309288
,1.86479609019
,1.92289164646
,1.98125692852
,2.03989809587
,2.09882044138
,2.1580285091
,2.21752619349
,2.2773168235
,2.33740323411
,2.39778782755
,2.45847262587
,2.51945931629
,2.58074929044
,2.6423436785
,2.70424337889
,2.76644908432
,2.82896130457
,2.89178038653
,2.95490653189
,3.01833981278
,3.08208018557
,3.14612750322
,3.21048152616
,3.27514193206
,3.3401083245
,3.40538024081
,3.47095715898
,3.53683850394
,3.60302365319
,3.66951194185
,3.73630266721
,3.80339509285
,3.87078845233
,3.93848195259
,4.00647477693
,4.0747660878
,4.14335502927
,4.21224072931
,4.28142230181
,4.35089884848
,4.42066946052
,4.49073322016
,4.56108920208
,4.63173647466
,4.70267410118
,4.77390114082
,4.84541664966
,4.91721968157
,4.98930928897
,5.06168452358
,5.13434443707
,5.20728808169
,5.28051451077
,5.35402277929
,5.42781194424
,5.50188106512
,5.57622920427
,5.65085542721
,5.72575880295
,5.8009384043
,5.87639330806
,5.9521225953
,6.02812535153
,6.10440066692
,6.1809476364
,6.25776535989
,6.33485294237
,6.412209494
,6.48983413025
,6.56772597196
,6.64588414544
,6.72430778253
,6.80299602064
,6.88194800282
,6.96116287779
,7.04063979998
,7.12037792952
,7.20037643233
,7.28063448005
,7.36115125011
,7.4419259257]

pl.figure()
pl.plot(A_dep,totals)
pl.show()












































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


