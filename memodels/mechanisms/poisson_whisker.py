import numpy as np
import math
from NeuroTools import stgen

def createRate(d, t):
    rate = []
    rater = []
    rate_sub2 = []
    rater_sub2 = []
    raterr_suball = []
    raterrr_suball = []
    for t in range(0,t):
        ts = 0
        tau1 = 10
        tau2 = 22
        alpha = 10.5
        beta = 10
        PrV1mean = math.pi
        PrV2mean = 5*(math.pi/4)
        PrV3mean = 0
        PrV4mean = (math.pi/4)
        PrV5mean = (math.pi/2)
        PrV6mean = 3*(math.pi/2)
        std = 69.53*(math.pi/180)

        G = 156.49485724*math.exp(-0.5*((d-PrV1mean)/std)**2) #156.49485724
        Gr = 156.49485724*math.exp(-0.5*((d-PrV2mean)/std)**2) #96.6017644
        G_sub2 = 156.49485724*math.exp(-0.5*((d-PrV3mean)/std)**2) #156.49485724
        Gr_sub2 = 156.49485724*math.exp(-0.5*((d-PrV4mean)/std)**2) #96.6017644

        Grr_suball = 156.49485724*math.exp(-0.5*((d-PrV5mean)/std)**2) #156.49485724
        Grrr_suball = 156.49485724*math.exp(-0.5*((d-PrV6mean)/std)**2) #96.6017644

        if (t-ts)<=0:
            thet = 0
        elif (t-ts)>0:
            thet = 1
            
        F = thet*alpha*(1-math.exp(-(t-ts)/tau1))-beta*(1-math.exp(-(t-ts)/tau2))
        R = G*F
        R_sub2 = G_sub2*F
        Rr = Gr*F
        Rr_sub2 = Gr_sub2*F

        Rrr_suball = Grr_suball*F
        Rrrr_suball = Grrr_suball*F

        # R = G
        # R_sub2 = G_sub2
        # Rr = Gr
        # Rr_sub2 = Gr_sub2
        rate.append(R)
        rater.append(Rr)
        rate_sub2.append(R_sub2)
        rater_sub2.append(Rr_sub2)

        raterr_suball.append(Rrr_suball)
        raterrr_suball.append(Rrrr_suball)

    rate_array = np.asarray(rate)
    rate_arrayr = np.asarray(rater)
    rate_array_sub2 = np.asarray(rate_sub2)
    rate_arrayr_sub2 = np.asarray(rater_sub2)

    rate_arrayrr_suball = np.asarray(raterr_suball)
    rate_arrayrrr_suball = np.asarray(raterrr_suball)

    return rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2, rate_arrayrr_suball, rate_arrayrrr_suball

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
    subnetworkall_PSarrayrr = []
    subnetworkall_PSarrayrrr = []
    for d in np.arange(-math.pi,2*math.pi,math.pi/16):
        # first make subnetwork spot inputs. 200ms delay, first 400ms is full field and then next 400ms is the spot of diameter d
        t = 400
        rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2, rate_arrayrr_suball, rate_arrayrrr_suball = createRate(d, t)
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
        if np.any(rate_arrayrr_suball==-0):
            for i in range(t):
                if rate_arrayrr_suball[i]==-0:
                    rate_arrayrr_suball[i]=rate_arrayrr_suball[i]*-1
        if np.any(rate_arrayrrr_suball==-0):
            for i in range(t):
                if rate_arrayrrr_suball[i]==-0:
                    rate_arrayrrr_suball[i]=rate_arrayrrr_suball[i]*-1       

        PS = createPoissonInputs(rate_array, t)
        PSr = createPoissonInputs(rate_arrayr, t)
        PS_sub2 = createPoissonInputs(rate_array_sub2,t)
        PSr_sub2 = createPoissonInputs(rate_arrayr_sub2,t)
        PSrr_suball = createPoissonInputs(rate_arrayrr_suball,t)
        PSrrr_suball = createPoissonInputs(rate_arrayrrr_suball,t)
        PS_delay = PS + 600
        PS_delayr = PSr + 600
        PS_delay_sub2 = PS_sub2 + 600
        PS_delayr_sub2 = PSr_sub2 + 600
        PS_delayrr_suball = PSrr_suball + 600
        PS_delayrrr_suball = PSrrr_suball + 600
    
        t = 400
        d = math.pi/2 # low stimulus
        rate_array,rate_arrayr,rate_array_sub2,rate_arrayr_sub2, rate_arrayrr_suball, rate_arrayrrr_suball = createRate(d,t)
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
        if np.any(rate_arrayrr_suball==-0):
            for i in range(t):
                if rate_arrayrr_suball[i]==-0:
                    rate_arrayrr_suball[i]=rate_arrayrr_suball[i]*-1
        if np.any(rate_arrayrrr_suball==-0):
            for i in range(t):
                if rate_arrayrrr_suball[i]==-0:
                    rate_arrayrrr_suball[i]=rate_arrayrrr_suball[i]*-1                         
        PS = createPoissonInputs(rate_array, t)
        PSr = createPoissonInputs(rate_arrayr, t)
        PS_sub2 = createPoissonInputs(rate_array_sub2,t)
        PSr_sub2 = createPoissonInputs(rate_arrayr_sub2,t)
        PSrr_suball = createPoissonInputs(rate_arrayrr_suball,t)
        PSrrr_suball = createPoissonInputs(rate_arrayrrr_suball,t)
        fullfield_delay = PS + 200
        fullfield_delayr = PSr + 200
        fullfield_delay_sub2 = PS_sub2 + 200
        fullfield_delayr_sub2 = PSr_sub2 + 200
        fullfield_delayrr_suball = PSrr_suball + 200
        fullfield_delayrrr_suball = PSrrr_suball + 200
    
        subnetwork1_PS = np.hstack((fullfield_delay, PS_delay))
        subnetwork1_PSr = np.hstack((fullfield_delayr, PS_delayr))
        subnetwork2_PS = np.hstack((fullfield_delay_sub2, PS_delay_sub2))
        subnetwork2_PSr = np.hstack((fullfield_delayr_sub2, PS_delayr_sub2))
        subnetworkall_PSrr = np.hstack((fullfield_delayrr_suball, PS_delayrr_suball))
        subnetworkall_PSrrr = np.hstack((fullfield_delayrrr_suball, PS_delayrrr_suball))
    
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
        subnetworkall_PSarrayrr.append(subnetworkall_PSrr)
        subnetworkall_PSarrayrrr.append(subnetworkall_PSrrr)

    subnetwork1_input = np.array(subnetwork1_PSarray)
    subnetwork2_input = np.array(subnetwork2_PSarray)
    subnetwork1_inputr = np.array(subnetwork1_PSarrayr)
    subnetwork2_inputr = np.array(subnetwork2_PSarrayr)
    subnetworkall_inputrr = np.array(subnetworkall_PSarrayrr)
    subnetworkall_inputrrr = np.array(subnetworkall_PSarrayrrr)

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
    np.save('../../spike_input/' + parent_dir_str+'/'+'extra_input1/subnetworkall_input' + '_' + trial_str + '.npy', subnetworkall_inputrr, allow_pickle = True)
    np.save('../../spike_input/' + parent_dir_str+'/'+'extra_input2/subnetworkall_input' + '_' + trial_str + '.npy', subnetworkall_inputrrr, allow_pickle = True)
    return subnetwork1_input, subnetwork2_input, subnetwork1_inputr, subnetwork2_inputr

