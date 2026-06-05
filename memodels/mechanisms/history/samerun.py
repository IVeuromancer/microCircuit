from neuron import h
from neuron.units import ms
from subnetwork import Subnetwork
import numpy as np
import time
import datetime
# import winsound

start = time.time()
starttime = datetime.datetime.now().time()
print ("start time :", starttime)

new_subnetwork1 = np.load('../../spike_input/subnetwork1_input.npy', allow_pickle = True)
new_subnetwork2 = np.load('../../spike_input/subnetwork2_input.npy', allow_pickle = True)
# print(new_subnetwork1.shape)
# print(new_subnetwork2.shape)
# print(new_subnetwork1) #print just to confirm it's the same as what was created in createSubnetworkInputs()
# print(new_subnetwork2) #print just to confirm it's the same as what was created in createSubnetworkInputs()

N=2
subnetwork = Subnetwork(0,N,100)
subnetwork2 = Subnetwork(5,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)

sameIN_time_series = []
sameIN_spike_events = []
for x in range(0,len(new_subnetwork1)): #loop through different spot diameter spike trains
    # change spike inputs to network
    subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
    subnetwork._vecstim.play(subnetwork.train_vec)
    subnetwork2.train_vec = h.Vector(new_subnetwork2[x,])
    subnetwork2._vecstim.play(subnetwork2.train_vec)
    # here is where I connect the two subnetworks with each other
    # connect two subnetworks by adding IN-IN axons and dendrites
    syn_IN1axon_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    syn_IN1axon_IN2.tau1 = 0.7
    syn_IN1axon_IN2.tau2 = 4.2
    syn_IN1axon_IN2.e = -80
    nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
    nc_IN1axon_IN2.weight[0] = 0.006
    nc_IN1axon_IN2.delay = 1

    syn_IN2axon_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    syn_IN2axon_IN1.tau1 = 0.7
    syn_IN2axon_IN1.tau2 = 4.2
    syn_IN2axon_IN1.e = -80
    nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
    nc_IN2axon_IN1.weight[0] = 0.006
    nc_IN2axon_IN1.delay = 1

    syn_IN1dend_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    syn_IN1dend_IN2.tau1 = 0.7
    syn_IN1dend_IN2.tau2 = 4.2
    syn_IN1dend_IN2.e = -80
    nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[17](0.9)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[17])
    nc_IN1dend_IN2.weight[0] = 0.006
    nc_IN1dend_IN2.delay = 1
    nc_IN1dend_IN2.threshold = -15

    syn_IN2dend_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    syn_IN2dend_IN1.tau1 = 0.7
    syn_IN2dend_IN1.tau2 = 4.2
    syn_IN2dend_IN1.e = -80
    nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[17](0.9)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[17])
    nc_IN2dend_IN1.weight[0] = 0.006
    nc_IN2dend_IN1.delay = 1
    nc_IN2dend_IN1.threshold = -15  
    
    # add Br inputs to own IN (sameIN network arrangement)
    synBr1toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    synBr1toIN1.tau1 = 1.6 * ms
    synBr1toIN1.tau2 = 3.6 * ms
    synBr1toIN1.e = 10 
    nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
    nc_Br1toIN1.delay = 1 * ms
    nc_Br1toIN1.weight[0] = 0.008

    synBr1toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[17](0.9))
    synBr1toIN2.tau1 = 0.3 * ms
    synBr1toIN2.tau2 = 2 * ms
    synBr1toIN2.e = 10
    nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
    nc_Br1toIN2.delay = 1 * ms
    nc_Br1toIN2.weight[0] = 0.008

    synBr2toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    synBr2toIN2.tau1 = 1.6 * ms
    synBr2toIN2.tau2 = 3.6 * ms
    synBr2toIN2.e = 10
    nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
    nc_Br2toIN2.delay = 1 * ms
    nc_Br2toIN2.weight[0] = 0.008

    synBr2toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[17](0.9))
    synBr2toIN1.tau1 = 0.3 * ms
    synBr2toIN1.tau2 = 2 * ms
    synBr2toIN1.e = 10
    nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
    nc_Br2toIN1.delay = 1 * ms
    nc_Br2toIN1.weight[0] = 0.008  
    
    t = h.Vector().record(h._ref_t)
    subnetwork._netstim.number = 0
    subnetwork2._netstim.number = 0
    h.v_init = -64
    h.celsius = 34
    h.tstop = 1000
    h.init()
    h.run()
    samevertstack = []
    samevertstack = np.vstack((
    subnetwork.IN_1.soma_v,
    subnetwork.IN_1.dend_prox_v,
    subnetwork.IN_1.dend_dist_v,
    subnetwork.IN_1.dend_dist2_v,
    subnetwork.TCcells[0].soma_v,
    subnetwork2.IN_1.soma_v,
    subnetwork2.IN_1.dend_prox_v,
    subnetwork2.IN_1.dend_dist_v,
    subnetwork2.IN_1.dend_dist2_v,
    subnetwork2.TCcells[0].soma_v))
    samevertstack2 = []
    samevertstack2 = [
    np.array(subnetwork.IN_1.spike_times),
    np.array(subnetwork.IN_1.release_times),
    np.array(subnetwork.TCcells[0].spike_times),
    np.array(subnetwork2.IN_1.spike_times),
    np.array(subnetwork2.IN_1.release_times),
    np.array(subnetwork2.TCcells[0].spike_times)]
    
    sameIN_time_series.append(samevertstack)
    sameIN_spike_events.append(samevertstack2)

sameIN_time_series = np.array(sameIN_time_series)
sameIN_spike_events = np.array(sameIN_spike_events)

print(sameIN_time_series.shape)
print(sameIN_spike_events.shape)

# Write the arrays to disk
datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
filename_timeseries = 'sameIN_time_series'
filename_spikes = 'sameIN_spike_events' 
txt_ext = '.txt'
pickle_ext = '.npy'
folder = '../../spike_output/'

with open(folder+filename_timeseries+datetimestr+txt_ext, 'w') as outfile:
    outfile.write('# Array shape: {0}\n'.format(sameIN_time_series.shape))
    for data_slice in sameIN_time_series:
        np.savetxt(outfile, data_slice, fmt='%-12.8f')
        outfile.write('# New diameter\n')   
        
np.save(folder+filename_spikes+datetimestr+pickle_ext, sameIN_spike_events, allow_pickle = True)

# duration = 5000  # milliseconds
# freq = 440  # Hz
# winsound.Beep(freq, duration)  
 
print(folder+filename_timeseries+datetimestr+txt_ext)
print(folder+filename_spikes+datetimestr+pickle_ext)

endtime = datetime.datetime.now().time()
print ("end time :", endtime)
end = time.time()
seconds_input = end - start
conversion = datetime.timedelta(seconds=seconds_input)
converted_time = str(conversion)
print(seconds_input)
print(converted_time)




