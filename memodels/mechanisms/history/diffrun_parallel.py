from neuron import h
from neuron.units import ms, mV
from subnetwork_parallel import Subnetwork
import numpy as np
import time
import datetime
import winsound

start = time.time()
starttime = datetime.datetime.now().time()
print ("start time :", starttime)
pc = h.ParallelContext()

new_subnetwork1 = np.load('../../spike_input/subnetwork1_input.npy', allow_pickle = True)
new_subnetwork2 = np.load('../../spike_input/subnetwork2_input.npy', allow_pickle = True)
# print(new_subnetwork1.shape)
# print(new_subnetwork2.shape)
# print(new_subnetwork1) #print just to confirm it's the same as what was created in createSubnetworkInputs()
# print(new_subnetwork2) #print just to confirm it's the same as what was created in createSubnetworkInputs()

pc.gid_clear() 
N = 2
subnetwork = Subnetwork(0,N,100)
subnetwork2 = Subnetwork(5,N,100)

diffIN_time_series = []
diffIN_spike_events = []
for x in range(0,len(new_subnetwork1)): #loop through different spot diameter spike trains
    # change spike inputs to network
    subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
    subnetwork._vecstim.play(subnetwork.train_vec)
    subnetwork2.train_vec = h.Vector(new_subnetwork2[x,])
    subnetwork2._vecstim.play(subnetwork2.train_vec)
    
    # connect two subnetworks by adding IN-IN axons and dendrites
    syn_IN1axon_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    syn_IN1axon_IN2.tau1 = 0.7
    syn_IN1axon_IN2.tau2 = 4.2
    syn_IN1axon_IN2.e = -80
    nc_IN1axon_IN2 = pc.gid_connect(subnetwork.IN_1._gid, syn_IN1axon_IN2)
    nc_IN1axon_IN2.weight[0] = 0.006
    nc_IN1axon_IN2.delay = 1
    nc_IN1axon_IN2.threshold = -15

    syn_IN2axon_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    syn_IN2axon_IN1.tau1 = 0.7
    syn_IN2axon_IN1.tau2 = 4.2
    syn_IN2axon_IN1.e = -80
    nc_IN2axon_IN1 = pc.gid_connect(subnetwork2.IN_1._gid, syn_IN2axon_IN1)
    nc_IN2axon_IN1.weight[0] = 0.006
    nc_IN2axon_IN1.delay = 1
    nc_IN2axon_IN1.threshold = -15

    syn_IN1dend_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    syn_IN1dend_IN2.tau1 = 0.7
    syn_IN1dend_IN2.tau2 = 4.2
    syn_IN1dend_IN2.e = -80
    nc_IN1dend_IN2 = pc.gid_connect(subnetwork.IN_1._gid3, syn_IN1dend_IN2)
    nc_IN1dend_IN2.weight[0] = 0.006
    nc_IN1dend_IN2.delay = 1
    nc_IN1dend_IN2.threshold = -15

    syn_IN2dend_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    syn_IN2dend_IN1.tau1 = 0.7
    syn_IN2dend_IN1.tau2 = 4.2
    syn_IN2dend_IN1.e = -80
    nc_IN2dend_IN1 = pc.gid_connect(subnetwork2.IN_1._gid3, syn_IN2dend_IN1)
    nc_IN2dend_IN1.weight[0] = 0.006
    nc_IN2dend_IN1.delay = 1
    nc_IN2dend_IN1.threshold = -15  
    
    # add Br inputs to other INs (diffIN network arrangement)
    synBr1toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.1))
    synBr1toIN2.tau1 = 1.6 * ms
    synBr1toIN2.tau2 = 3.6 * ms
    synBr1toIN2.e = 10
    nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
    nc_Br1toIN2.delay = 1 * ms
    nc_Br1toIN2.weight[0] = 0.008

    synBr1toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[17](0.9))
    synBr1toIN1.tau1 = 0.3 * ms
    synBr1toIN1.tau2 = 2 * ms
    synBr1toIN1.e = 10
    nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
    nc_Br1toIN1.delay = 1 * ms
    nc_Br1toIN1.weight[0] = 0.008

    synBr2toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.1))
    synBr2toIN1.tau1 = 1.6 * ms
    synBr2toIN1.tau2 = 3.6 * ms
    synBr2toIN1.e = 10
    nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
    nc_Br2toIN1.delay = 1 * ms
    nc_Br2toIN1.weight[0] = 0.008

    synBr2toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[17](0.9))
    synBr2toIN2.tau1 = 0.3 * ms
    synBr2toIN2.tau2 = 2 * ms
    synBr2toIN2.e = 10
    nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
    nc_Br2toIN2.delay = 1 * ms
    nc_Br2toIN2.weight[0] = 0.008
    
    # the TCcells don't synapse onto anything so they don't need to send spikes to all the other processors
    # pc.outputcell(0) 
    # pc.outputcell(1) 
    # pc.outputcell(2) 
    # pc.outputcell(5) 
    # pc.outputcell(6) 
    # pc.outputcell(7) 
    
    t = h.Vector().record(h._ref_t)
    subnetwork._netstim.number = 0 #remove netstim because we just want vecstim from the poisson inputs we created
    subnetwork2._netstim.number = 0
    pc.set_maxstep(10 * ms)
    h.v_init = -64
    h.celsius = 34
    h.init()
    pc.psolve(1000 * ms)
    
    local_data = {cell._gid: list(cell.soma_v) for cell in subnetwork.TCcells}
    # local_data.update = {cell._gid: list(cell.soma_v) for cell in subnetwork2.TCcells}
    all_data = pc.py_alltoall([local_data] + [None] * (pc.nhost() - 1))

    if pc.id() == 0:
        # combine the data from the various processes
        data = {}
        for process_data in all_data:
            data.update(process_data)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        filename_timeseries = 'diffIN_time_series'
        # filename_spikes = 'diffIN_spike_events' 
        txt_ext = '.txt'
        pickle_ext = '.npy'
        folder = '../../spike_output/'
        np.save(folder+filename_timeseries+datetimestr+pickle_ext, data, allow_pickle = True)
     
        print(folder+filename_timeseries+datetimestr+txt_ext)
        # print(folder+filename_spikes+datetimestr+pickle_ext)
        
frequency = 440
duration = 5000
winsound.Beep(frequency, duration)

endtime = datetime.datetime.now().time()
print ("end time :", endtime)
end = time.time()
seconds_input = end - start
conversion = datetime.timedelta(seconds=seconds_input)
converted_time = str(conversion)
print(seconds_input)
print(converted_time)


pc.barrier()
pc.done()
h.quit()



