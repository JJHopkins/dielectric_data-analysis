#!/opt/local/bin/python

# Import commands
import matplotlib.pyplot as pl
from numpy import *
from scipy.optimize import leastsq

def decaying_sine(amplitude,frequency,tau,t0,times):
    """
    Decaying sinusoid
    """
    signal = amplitude*sin(2*pi*frequency*(times-t0))*exp(-(times-t0)/tau)
    signal[times<t0]=0

    return(signal)

def residuals(params,y,x):
    a,f,tau,t0 = params
    delta = y-decaying_sine(a,f,tau,t0,x)
    return delta


# Generate time-axis
num_points = 1e3
time_stamps = linspace(0,10,num_points)

# signal parameters
amp,freq,tau,t0 = 1.73,2.69,2,2.5
true_signal = decaying_sine(amp,freq,tau,t0,time_stamps)

# Add some noise
noise = 10*random.normal(0,0.05,len(true_signal))
measured_signal = true_signal+noise

# initial guess
params0 = [1,1,1,2.5]

params_estimate,success = leastsq(residuals, params0,
                          args=(measured_signal,time_stamps))


print 'simulated: A=%1.5f f=%1.5f tau=%1.5f t0=%1.5f'%(
    amp,freq,tau,2.5)
print 'recovered: A=%1.5f f=%1.5f tau=%1.5f t0=%1.5f'%(
    params_estimate[0],params_estimate[1],params_estimate[2],
    params_estimate[3])

recovered_signal=decaying_sine(params_estimate[0],
                               params_estimate[1],
                               params_estimate[2],
                               params_estimate[3],
                               time_stamps)


# Make figure
pl.figure()
#pl.subplot(121)
pl.plot(time_stamps,measured_signal,label='data')
pl.plot(time_stamps,true_signal,'r',label='true')
pl.plot(time_stamps,recovered_signal,label='recovered')
pl.legend()
#pl.subplot(122)
#pl.plot(time_stamps,true_signal-recovered_signal,'g')
#pl.title('residuals')
pl.show()





