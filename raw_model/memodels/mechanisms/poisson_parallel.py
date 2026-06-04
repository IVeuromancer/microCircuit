import numpy as np
import math
from NeuroTools import stgen

def createRate(d, t):
    rate = []
    for t in range(0,t):
        w = 0.85
        #lbkg(1-w) = 36.8
        lbkg = 36.8/0.15
        #lspot(1-w) = 56.5
        lspot = 56.5/0.15
        a1 = 0.62
        a2 = 1.26
        tau1 = 10
        tau2 = 22
        alpha = 12
        beta = 11.26
        ts = 0
        
        # Gprime = lbkg(1-w)+(lspot-lbkg)(1-math.exp(-d**2/(4*a1**2))-w(1-math.exp(-d**2/(4*a2**2))))
        Gprime = 36.8 + (lspot-lbkg)*(1-math.exp(-d**2/(4*a1**2))-w*(1-math.exp(-d**2/(4*a2**2))))

        if Gprime<0:
            thet = 0
        elif Gprime>0:
            thet = 1

        G = Gprime*thet

        if (t-ts)<=0:
            thet2 = 0
        elif (t-ts)>0:
            thet2 = 1

        F = thet2*alpha*(1-math.exp(-(t-ts)/tau1))-beta*(1-math.exp(-(t-ts)/tau2))
        R = G*F
        rate.append(R)
    rate_array = np.asarray(rate)
    return rate_array
    
def createPoissonInputs(rate_array, t):
    # NOTE: HAD TO UPDATE STGEN.PY FILE BECAUSE LATEST VERSION OF NEUROTOOLS IS NOT COMPATIBLE WITH LATEST VERISON OF PYTHON/NUMPY
    # specifically, I had to change xrange() to range() and add () to print. Also update relative import paths to the signals folder
    # I also had to change remove all float numbers to integers and add .astype(int) to a few lines for numpy compatibility. also updated 
    # random number generator
    st_gen = stgen.StGen()
    intervals=np.arange(0,t)
    t_stop = t
    PS = st_gen.inh_poisson_generator(rate_array,intervals, t_stop, array = True) # important to have rate in Hz and all other times in ms
    return PS

def createSubnetworkInputs():

    subnetwork1_PSarray = []
    subnetwork2_PSarray = []
    for d in np.arange(1,2.5,0.5):
        # first make subnetwork spot inputs. 200ms delay, first 400ms is full field and then next 400ms is the spot of diameter d
        t = 400
        rate_array = createRate(d, t)
        PS = createPoissonInputs(rate_array, t)
        PS_delay = PS + 600
    
        t = 400
        d = 0 # full field
        rate_array = createRate(d,t)
        PS = createPoissonInputs(rate_array, t)
        fullfield_delay = PS + 200
    
        subnetwork1_PS = np.hstack((fullfield_delay, PS_delay))
    
        # make subnetwork2 full field for 800ms after 200ms delay
        t = 800
        d = 0
        rate_array = createRate(d, t)
        PS = createPoissonInputs(rate_array, t)
        subnetwork2_PS = PS + 200
        
        subnetwork1_PSarray.append(subnetwork1_PS)
        subnetwork2_PSarray.append(subnetwork2_PS)
        
    subnetwork1_input = np.array(subnetwork1_PSarray)
    subnetwork2_input = np.array(subnetwork2_PSarray)
    print(subnetwork1_input.shape)
    print(subnetwork2_input.shape)
    print(subnetwork1_input)
    print(subnetwork2_input)
    
    np.save('../../spike_input/subnetwork1_input.npy', subnetwork1_input, allow_pickle = True)
    np.save('../../spike_input/subnetwork2_input.npy', subnetwork2_input, allow_pickle = True)
    return subnetwork1_input, subnetwork2_input

createSubnetworkInputs()    