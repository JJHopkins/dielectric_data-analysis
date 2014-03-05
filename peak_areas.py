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

import matplotlib               
import numpy                    # http://numpy.scipy.org/
import scipy                    # http://scipy.org/
import scipy.optimize, scipy.special, scipy.stats
from matplotlib import pyplot
import scipy.integrate as grt

def lorentzian(x,*p) :
    # A lorentzian peak with:
    #   Peak height above background : p[0]
    #   Central value                : p[1]
    #   Full Width at Half Maximum   : p[2]
    #return (p[0]/numpy.pi)/(1.0+((x-p[1])/p[2])**2)
    return (p[0])/(1.0+((x-p[1])/p[2])**2)

#---------------------------------------------------------
#
# Fit MEASURED peaks
#Y24903L_pk.PRN
x_data_expt, y_data_expt = numpy.loadtxt('data/Y24903L.PRN', unpack=True, usecols = [0,1])
#x_data_expt, y_data_expt = numpy.loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
#x_data_expt, y_data_expt = numpy.loadtxt('expt_viacsv.dat', unpack=True, usecols = [0,1])

# 1st peak
xMin0 = 0.0
xMax0 = 10.5
p0_guess = (3.91360,10.3022,0.40694)

y_fitData=y_data_expt[(x_data_expt>xMin0)*(x_data_expt<xMax0)]
x_fitData=x_data_expt[(x_data_expt>xMin0)*(x_data_expt<xMax0)]

exptPeak0Params, exptPeak0CovMat = scipy.optimize.curve_fit(lorentzian,
        x_fitData, y_fitData, p0=p0_guess)

exptPeak0Lorentz=lorentzian(x_data_expt,*exptPeak0Params)

# 2nd peak
xMin1 = 11.0
xMax1 = 13.0
p1_guess = (2.53875,11.7171,1.45631)

y_fitData=y_data_expt[(x_data_expt>xMin1)*(x_data_expt<xMax1)]
x_fitData=x_data_expt[(x_data_expt>xMin1)*(x_data_expt<xMax1)]

exptPeak1Params, exptPeak1CovMat = scipy.optimize.curve_fit(lorentzian,
        x_fitData, y_fitData, p0=p1_guess)

exptPeak1Lorentz=lorentzian(x_data_expt,*exptPeak1Params)

# 3rd peak
xMin2 = 13.0
xMax2 = 14.5
p2_guess = (2.22888,13.7171,1.14198)

y_fitData=y_data_expt[(x_data_expt>xMin2)*(x_data_expt<xMax2)]
x_fitData=x_data_expt[(x_data_expt>xMin2)*(x_data_expt<xMax2)]

exptPeak2Params, exptPeak2CovMat = scipy.optimize.curve_fit(lorentzian,
        x_fitData, y_fitData, p0=p1_guess)

exptPeak2Lorentz=lorentzian(x_data_expt,*exptPeak2Params)

#---------------------------------------------------------
#
# Fit CALCULATED data
#
#x_data_theory, y_data_theory = numpy.loadtxt('iwyc_a_input.dat', unpack=True, usecols = [0,1])
x_data_theory, y_data_theory = numpy.loadtxt('data/Y23605L.PRN', unpack=True, usecols = [0,1])

# NOTE: the theory data does not /have/ a zeroth peak to match the experimental
# data!

# 2nd peak
xMin1 = 8.5
xMax1 = 10.5
p1_guess = (2.53875,11.7171,1.45631)

y_fitData=y_data_theory[(x_data_theory>xMin1)*(x_data_theory<xMax1)]
x_fitData=x_data_theory[(x_data_theory>xMin1)*(x_data_theory<xMax1)]

theoryPeak1Params, theoryPeak1CovMat = scipy.optimize.curve_fit(lorentzian,
        x_fitData, y_fitData, p0=p1_guess)

# Now, shift (translate) the theoretical data such that the 2nd peak has the same location
# as the experimental data
abscissaShift = exptPeak1Params[1] - theoryPeak1Params[1]
x_data_theory += abscissaShift

# don't forget that we've shifted the center of the lorentzian fit!
theoryPeak1Params[1] += abscissaShift

theoryPeak1Lorentz=lorentzian(x_data_theory,*theoryPeak1Params)

# Renormalise theoretical data so that the area under Lorentz peak 1 is the same
# in theory and experiment:
exptPeak1LorentzArea=numpy.trapz(exptPeak1Lorentz,x_data_expt)
theoryPeak1LorentzArea=numpy.trapz(theoryPeak1Lorentz,x_data_theory)

# XXX: actually, using the area is pretty bad
#y_data_theory *= exptPeak1LorentzArea/theoryPeak1LorentzArea
#theoryPeak1Params[0] *= exptPeak1LorentzArea/theoryPeak1LorentzArea

# Just take the ratios of the peak heights
y_data_theory *= max(exptPeak1Lorentz)/max(theoryPeak1Lorentz)
theoryPeak1Params[0] *= max(exptPeak1Lorentz)/max(theoryPeak1Lorentz)

# recompute the fit for plotting
theoryPeak1Lorentz=lorentzian(x_data_theory,*theoryPeak1Params)

#----------------------------------------------------------------------------
#
# Synthesise a new theoretical spectrum with the addition of the exciton peak
#

# make a copy of the theoretical data (copying is a good idea; setting things
# equal to each other is extremely dangerous!)
y_data_synth = numpy.copy(y_data_theory)
x_data_synth = x_data_theory

# Re-compute expt Lorentz peak0 using theoretical abscissa
exptPeak0LorentzSynth=lorentzian(x_data_theory,*exptPeak0Params)

# Locate point of intersection between expt Lorentz peak 1 and calculated data
xSwitchover = x_data_theory[exptPeak0LorentzSynth>y_data_theory][-1]

# replace entries of theory calculation < intersection point with the fitted
# Lorentzian peak 0
y_data_synth[x_data_synth<=xSwitchover] = \
    exptPeak0LorentzSynth[x_data_synth<=xSwitchover]


#---------------------------------------------------------
#
# PLOT!!!
#


# Plot the experimental data
pyplot.figure()

# Plot the lorentzian fits to the peaks in the experimental data
pyplot.plot(x_data_expt,exptPeak0Lorentz, color = 'r', 
    label='Lorentzian Fit (expt peak 0)', linestyle='--')
pyplot.plot(x_data_expt,exptPeak1Lorentz, color = 'r', 
    label='Lorentzian Fit (expt peak 1)', linestyle='--')
pyplot.plot(x_data_expt,exptPeak2Lorentz, color = 'r', 
    label='Lorentzian Fit (expt peak 2)', linestyle='--')

pyplot.plot(x_data_expt,y_data_expt, color = 'k', linestyle = 'none',
        label='measured', marker=',')

pyplot.legend()
pyplot.title('Experimental Data & Fits')
pyplot.show()

# Plot the theory/calculated data
pyplot.figure()
pyplot.plot(x_data_theory, y_data_theory, color = 'k', label='theory')

# Plot the lorentzian fits to the peaks in the experimental data
pyplot.plot(x_data_theory,theoryPeak1Lorentz, color = 'g', 
    label='Lorentzian Fit (theory peak 1)', linestyle='--')

pyplot.title('Theoretical Data & Fits')

pyplot.legend()
pyplot.show()

# Finally, plot the newly synthesised spectrum next to the original data for
# visual inspection...
pyplot.figure()
pyplot.plot(x_data_synth,y_data_synth, color='r', label='synthesised')
pyplot.plot(x_data_expt,y_data_expt, color = 'k', linestyle = 'none',
        label='measured', marker=',')
pyplot.axvline(xSwitchover,color='m',label='end of synthesis')
pyplot.legend()
pyplot.show()


