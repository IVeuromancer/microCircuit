import numpy as np
import math
from NeuroTools import stgen
import scipy.special


def createRate(d, t):
    rate = []
    rater = []
    rate_sub2 = []
    rater_sub2 = []
    for t in range(0,t):
        w = 0.85
        #lbkg(1-w) = 36.8
        lbkg = 36.8/0.15
        #lspot(1-w) = 56.5
        lspot = 56.5/0.15
        a1 = 0.62
        a2 = 1.26
        ra = 0.99
        tau1 = 10
        tau2 = 22
        alpha = 12
        beta = 11.26
        ts = 0
        
        # Gprime = lbkg(1-w)+(lspot-lbkg)(1-math.exp(-d**2/(4*a1**2))-w(1-math.exp(-d**2/(4*a2**2))))
        Gprime = 36.8 + (lspot-lbkg)*(1-math.exp(-(d**2)/(4*(a1**2)))-(0.85*(1-math.exp(-(d**2)/(4*(a2**2))))))
        Gprime_sub2 = 56.5 + (lbkg-lspot)*(1-math.exp(-(d**2)/(4*(a1**2)))-(0.85*(1-math.exp(-(d**2)/(4*(a2**2))))))

        if Gprime<=0:
            thet = 0
        elif Gprime>0:
            thet = 1

        G = Gprime*thet

        if Gprime_sub2<=0:
            thet2 = 0
        elif Gprime_sub2>0:
            thet2 = 1

        G_sub2 = Gprime_sub2*thet2

        def the_sum(d, a):
            return sum(
                (ra/a)**(2*m) / (math.factorial(m)) * scipy.special.gammainc(m + 1, (d**2)/(4*(a**2)))
                 for m in range(30)
            )

        sum1 = the_sum(d,a1)
        sum2 = the_sum(d,a2)
        Grprime = 36.8 + (lspot-lbkg)*(((math.exp(-(ra**2)/(a1**2)))*sum1)-(0.85*((math.exp(-(ra**2)/(a2**2))*sum2))))
        Grprime_sub2 = 56.5 + (lbkg-lspot)*(((math.exp(-(ra**2)/(a1**2)))*sum1)-(0.85*((math.exp(-(ra**2)/(a2**2))*sum2))))

        if Grprime<=0:
            thet3 = 0
        elif Grprime>0:
            thet3 = 1

        Gr = Grprime*thet3

        if Grprime_sub2<=0:
            thet4 = 0
        elif Grprime_sub2>0:
            thet4 = 1

        Gr_sub2 = Grprime_sub2*thet4

        if (t-ts)<=0:
            thet5 = 0
        elif (t-ts)>0:
            thet5 = 1

        F = thet5*alpha*(1-math.exp(-(t-ts)/tau1))-beta*(1-math.exp(-(t-ts)/tau2))
        R = G*F
        R_sub2 = G_sub2*F
        Rr = Gr*F
        Rr_sub2 = Gr_sub2*F
        rate.append(R)
        rater.append(Rr)
        rate_sub2.append(R_sub2)
        rater_sub2.append(Rr_sub2)
    rate_array = np.asarray(rate)
    rate_arrayr = np.asarray(rater)
    rate_array_sub2 = np.asarray(rate_sub2)
    rate_arrayr_sub2 = np.asarray(rater_sub2)
    return rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2

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

def createSubnetworkInputs(trial, parent_dir):

    subnetwork1_PSarray = []
    subnetwork2_PSarray = []
    subnetwork1_PSarrayr = []
    subnetwork2_PSarrayr = []
    for d in np.arange(0,10,0.1):
        # first make subnetwork spot inputs. 200ms delay, first 400ms is full field and then next 400ms is the spot of diameter d
        t = 400
        rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2 = createRate(d, t)
        if np.any(rate_array==-0):
            for i in range(t):
                if rate_array[i]==-0:
                    rate_array[i]=0.
        if np.any(rate_arrayr==-0):
            for i in range(t):
                if rate_arrayr[i]==-0:
                    rate_arrayr[i]=0.
        if np.any(rate_array_sub2==-0):
            for i in range(t):
                if rate_array_sub2[i]==-0:
                    rate_array_sub2[i]=rate_array_sub2[i]*-1
        if np.any(rate_arrayr_sub2==-0):
            for i in range(t):
                if rate_arrayr_sub2[i]==-0:
                    rate_arrayr_sub2[i]=rate_arrayr_sub2[i]*-1                           
        PS = createPoissonInputs(rate_array, t)
        PSr = createPoissonInputs(rate_arrayr, t)
        PS_sub2 = createPoissonInputs(rate_array_sub2,t)
        PSr_sub2 = createPoissonInputs(rate_arrayr_sub2,t)
        PS_delay = PS + 600
        PS_delayr = PSr + 600
        PS_delay_sub2 = PS_sub2 + 600
        PS_delayr_sub2 = PSr_sub2 + 600
    
        t = 400
        d = 10 # full field
        rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2 = createRate(d,t)
        if np.any(rate_array==-0):
            for i in range(t):
                if rate_array[i]==-0:
                    rate_array[i]=0.
        if np.any(rate_arrayr==-0):
            for i in range(t):
                if rate_arrayr[i]==-0:
                    rate_arrayr[i]=0.
        if np.any(rate_array_sub2==-0):
            for i in range(t):
                if rate_array_sub2[i]==-0:
                    rate_array_sub2[i]=rate_array_sub2[i]*-1
        if np.any(rate_arrayr_sub2==-0):
            for i in range(t):
                if rate_arrayr_sub2[i]==-0:
                    rate_arrayr_sub2[i]=rate_arrayr_sub2[i]*-1                      
        PS = createPoissonInputs(rate_array, t)
        PSr = createPoissonInputs(rate_arrayr, t)
        PS_sub2 = createPoissonInputs(rate_array_sub2,t)
        PSr_sub2 = createPoissonInputs(rate_arrayr_sub2,t)
        fullfield_delay = PS + 200
        fullfield_delayr = PSr + 200
        fullfield_delay_sub2 = PS_sub2 + 200
        fullfield_delayr_sub2 = PSr_sub2 + 200
    
        subnetwork1_PS = np.hstack((fullfield_delay, PS_delay))
        subnetwork1_PSr = np.hstack((fullfield_delayr, PS_delayr))
        subnetwork2_PS = np.hstack((fullfield_delay_sub2, PS_delay_sub2))
        subnetwork2_PSr = np.hstack((fullfield_delayr_sub2, PS_delayr_sub2))
    
        # make subnetwork2 full field for 800ms after 200ms delay
        # t = 800
        # d = 10
        # rate_array,rate_arrayr = createRate(d, t)
        # PS = createPoissonInputs(rate_array, t)
        # PSr = createPoissonInputs(rate_arrayr, t)
        # subnetwork2_PS = PS + 200
        # subnetwork2_PSr = PSr + 200
        
        subnetwork1_PSarray.append(subnetwork1_PS)
        subnetwork2_PSarray.append(subnetwork2_PS)
        subnetwork1_PSarrayr.append(subnetwork1_PSr)
        subnetwork2_PSarrayr.append(subnetwork2_PSr)

    subnetwork1_input = np.array(subnetwork1_PSarray)
    subnetwork2_input = np.array(subnetwork2_PSarray)
    subnetwork1_inputr = np.array(subnetwork1_PSarrayr)
    subnetwork2_inputr = np.array(subnetwork2_PSarrayr)

    # print(subnetwork1_input.shape)
    # print(subnetwork2_input.shape)
    # print(subnetwork1_input)
    # print(subnetwork2_input)
    trial_str = str(trial)
    parent_dir_str = str(parent_dir)
    np.save('../../spike_input/' + parent_dir_str+'/'+'Br1/subnetwork1_input' + '_' + trial_str + '.npy', subnetwork1_input, allow_pickle = True)
    np.save('../../spike_input/' + parent_dir_str+'/'+'Br1/subnetwork2_input' + '_' + trial_str + '.npy', subnetwork2_input, allow_pickle = True)
    np.save('../../spike_input/' + parent_dir_str+'/'+'Br2/subnetwork1_input' + '_' + trial_str + '.npy', subnetwork1_inputr, allow_pickle = True)
    np.save('../../spike_input/' + parent_dir_str+'/'+'Br2/subnetwork2_input' + '_' + trial_str + '.npy', subnetwork2_inputr, allow_pickle = True)
    return subnetwork1_input, subnetwork2_input, subnetwork1_inputr, subnetwork2_inputr

