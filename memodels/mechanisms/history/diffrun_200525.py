from neuron import h
from subnetwork_200525 import Subnetwork
import numpy as np
import time
import datetime
# import winsound

def runDiff(parent_dir, trial):

    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
    trial_str = str(trial)
    
    new_subnetwork1 = np.load('../../spike_input/subnetwork1_input' + '_' + trial_str + '.npy', allow_pickle = True)
    new_subnetwork2 = np.load('../../spike_input/subnetwork2_input' + '_' + trial_str + '.npy', allow_pickle = True)
    
    N=4
    subnetwork = Subnetwork(0,N,100)
    subnetwork2 = Subnetwork(9,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)
    
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): #loop through different spot diameter spike trains
        # change spike inputs to network
        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork2.train_vec = h.Vector(new_subnetwork2[x,])
        subnetwork2._vecstim.play(subnetwork2.train_vec)
        # here is where I connect the two subnetworks with each other
        # connect two subnetworks by adding IN-IN axons and dendrites
        syn_IN1axon_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1axon_IN2.tau1 = 0.71
        syn_IN1axon_IN2.tau2 = 4.18
        syn_IN1axon_IN2.e = -80
        nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
        nc_IN1axon_IN2.weight[0] = 0.0053
        nc_IN1axon_IN2.delay = 1
    
        syn_IN2axon_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2axon_IN1.tau1 = 0.71
        syn_IN2axon_IN1.tau2 = 4.18
        syn_IN2axon_IN1.e = -80
        nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
        nc_IN2axon_IN1.weight[0] = 0.0053
        nc_IN2axon_IN1.delay = 1
    
        syn_IN1dend_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1dend_IN2.tau1 = 0.71
        syn_IN1dend_IN2.tau2 = 4.18
        syn_IN1dend_IN2.e = -80
        nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[85](1)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[85])
        nc_IN1dend_IN2.weight[0] = 0.0053
        nc_IN1dend_IN2.delay = 0.1
        nc_IN1dend_IN2.threshold = -34
    
        syn_IN2dend_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2dend_IN1.tau1 = 0.71
        syn_IN2dend_IN1.tau2 = 4.18
        syn_IN2dend_IN1.e = -80
        nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[85](1)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[85])
        nc_IN2dend_IN1.weight[0] = 0.0053
        nc_IN2dend_IN1.delay = 0.1
        nc_IN2dend_IN1.threshold = -34
        
        # add Br inputs to other INs (diffIN network arrangement)
        synBr1toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        synBr1toIN2.tau1 = 0.32
        synBr1toIN2.tau2 = 2.08
        synBr1toIN2.e = 10
        nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        nc_Br1toIN2.delay = 1
        nc_Br1toIN2.weight[0] = 0.003
    
        synBr1toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[85](1))
        synBr1toIN1.tau1 = 0.29
        synBr1toIN1.tau2 = 2.15
        synBr1toIN1.e = 10
        nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        nc_Br1toIN1.delay = 1
        nc_Br1toIN1.weight[0] = 0.003
    
        synBr2toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        synBr2toIN1.tau1 = 0.32
        synBr2toIN1.tau2 = 2.08
        synBr2toIN1.e = 10
        nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        nc_Br2toIN1.delay = 1 
        nc_Br2toIN1.weight[0] = 0.003
    
        synBr2toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[85](1))
        synBr2toIN2.tau1 = 0.29
        synBr2toIN2.tau2 = 2.15
        synBr2toIN2.e = 10
        nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        nc_Br2toIN2.delay = 1 
        nc_Br2toIN2.weight[0] = 0.003
        
        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()
        #<3  diffvertstack = []
        #<3  diffvertstack = np.vstack((
        #<3  subnetwork.IN_1.soma_v,
        #<3  subnetwork.IN_1.dend_prox_v,
        #<3  subnetwork.IN_1.dend_dist_v,
        #<3  subnetwork.IN_1.dend_dist2_v,
        #<3  subnetwork.IN_1.dend_dist3_v,
        #<3  subnetwork.IN_1.dend_dist4_v,
        #<3  subnetwork.IN_1.dend_dist5_v,
        #<3  subnetwork.TCcells[0].soma_v,
        #<3  subnetwork.TCcells[1].soma_v,
        #<3  subnetwork.TCcells[2].soma_v,
        #<3  subnetwork.TCcells[3].soma_v,
        #<3  subnetwork2.IN_1.soma_v,
        #<3  subnetwork2.IN_1.dend_prox_v,
        #<3  subnetwork2.IN_1.dend_dist_v,
        #<3  subnetwork2.IN_1.dend_dist2_v,
        #<3  subnetwork2.IN_1.dend_dist3_v,
        #<3  subnetwork2.IN_1.dend_dist4_v,
        #<3  subnetwork2.IN_1.dend_dist5_v,
        #<3  subnetwork2.TCcells[0].soma_v,
        #<3  subnetwork2.TCcells[1].soma_v,
        #<3  subnetwork2.TCcells[2].soma_v,
        #<3  subnetwork2.TCcells[3].soma_v,))
        diffvertstack2 = []
        diffvertstack2 = [
        np.array(subnetwork.IN_1.spike_times),
        np.array(subnetwork.IN_1.release_times),
        np.array(subnetwork.IN_1.release_times2),
        np.array(subnetwork.IN_1.release_times3),
        np.array(subnetwork.IN_1.release_times4),
        np.array(subnetwork.IN_1.release_times5),
        np.array(subnetwork.TCcells[0].spike_times),
        np.array(subnetwork.TCcells[1].spike_times),
        np.array(subnetwork.TCcells[2].spike_times),
        np.array(subnetwork.TCcells[3].spike_times),
        np.array(subnetwork2.IN_1.spike_times),
        np.array(subnetwork2.IN_1.release_times),
        np.array(subnetwork2.IN_1.release_times2),
        np.array(subnetwork2.IN_1.release_times3),
        np.array(subnetwork2.IN_1.release_times4),
        np.array(subnetwork2.IN_1.release_times5),
        np.array(subnetwork2.TCcells[0].spike_times),
        np.array(subnetwork2.TCcells[1].spike_times),
        np.array(subnetwork2.TCcells[2].spike_times),
        np.array(subnetwork2.TCcells[3].spike_times)]
        
        #<3 diffIN_time_series.append(diffvertstack)
        #<3 diffIN_spike_events.append(diffvertstack2)
        #<3 diffIN_spike_events = np.array(diffIN_spike_events)
        diffIN_spike_events = np.array(diffvertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        #<3 filename_timeseries = 'diffIN_time_series'
        filename_spikes = 'diffIN_spike_events' 
        #<3 txt_ext = '.txt'
        pickle_ext = '.npy'
        folder = '../../spike_output/' + parent_dir + '/diffrun_' + trial_str + '/'
        xstr = str(x)
        np.save(folder+filename_spikes+datetimestr+'-'+xstr+pickle_ext, diffIN_spike_events, allow_pickle = True)

    #<3 diffIN_time_series = np.array(diffIN_time_series)
    
    #<3 print(diffIN_time_series.shape)
    #<3 print(diffIN_spike_events.shape)
    
    # Write the arrays to disk
    
    #<3 with open(folder+filename_timeseries+datetimestr+txt_ext, 'w') as outfile:
    #<3     outfile.write('# Array shape: {0}\n'.format(diffIN_time_series.shape))
    #<3     for data_slice in diffIN_time_series:
    #<3         np.savetxt(outfile, data_slice, fmt='%-12.8f')
    #<3         outfile.write('# New diameter\n')   
            
     
    #<3 print(folder+filename_timeseries+datetimestr+txt_ext)
    #<3 print(folder+filename_spikes+datetimestr+pickle_ext)

    
    endtime = datetime.datetime.now().time()
    print ("end time :", endtime)
    end = time.time()
    seconds_input = end - start
    conversion = datetime.timedelta(seconds=seconds_input)
    converted_time = str(conversion)
    # print(seconds_input)
    print(converted_time)
    
    # freq = 440
    # dur = 1000
    # winsound.Beep(freq,dur)




