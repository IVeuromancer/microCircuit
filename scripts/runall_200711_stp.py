from neuron import h
from subnetwork_200704_stp import Subnetwork
import poisson_whisker
import numpy as np
import time
import datetime
import winsound
import os
import multiprocessing
from multiprocessing import freeze_support, Pool
import convert_inputs


def diffrun(z):    
    parent_dir = 'run_200711_sameopp_stp'

    N=4
    subnetwork = Subnetwork(0,N,100)
    subnetwork2 = Subnetwork(9,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)

    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new_subnetwork2 = np.load('../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork2 = np.load('../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../spike_output/' + parent_dir + '/diffrun_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,16): #loop through different spot diameter spike trains
        # change spike inputs to network

        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        subnetwork2.train_vec = h.Vector(new_subnetwork2[x,])
        subnetwork2._vecstim.play(subnetwork2.train_vec)
        subnetwork2.train_vec2 = h.Vector(new2_subnetwork2[x,])
        subnetwork2._vecstim2.play(subnetwork2.train_vec2)

        # here is where I connect the two subnetworks with each other
        # connect two subnetworks by adding IN-IN axons and dendrites
        syn_IN1axon_IN2 = h.DetGABAAB(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1axon_IN2.tau_r_GABAA = 0.71
        syn_IN1axon_IN2.tau_d_GABAA = 4.18
        syn_IN1axon_IN2.e_GABAA = -80
        syn_IN1axon_IN2.Use = 0.25
        syn_IN1axon_IN2.Dep = 995
        syn_IN1axon_IN2.Fac = 3
        syn_IN1axon_IN2.u0 = 0.5
        nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
        nc_IN1axon_IN2.weight[0] = 10 #5
        nc_IN1axon_IN2.delay = 1
    
        syn_IN2axon_IN1 = h.DetGABAAB(subnetwork.IN_1.model.dend[25](1))
        syn_IN2axon_IN1.tau_r_GABAA = 0.71
        syn_IN2axon_IN1.tau_d_GABAA = 4.18
        syn_IN2axon_IN1.e_GABAA = -80
        syn_IN2axon_IN1.Use = 0.25
        syn_IN2axon_IN1.Dep = 995
        syn_IN2axon_IN1.Fac = 3
        syn_IN2axon_IN1.u0 = 0.5
        nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
        nc_IN2axon_IN1.weight[0] = 10 #5
        nc_IN2axon_IN1.delay = 1
    
        syn_IN1dend_IN2 = h.DetGABAAB(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1dend_IN2.tau_r_GABAA = 0.71
        syn_IN1dend_IN2.tau_d_GABAA = 4.18
        syn_IN1dend_IN2.e_GABAA = -80
        syn_IN1dend_IN2.Use = 0.25
        syn_IN1dend_IN2.Dep = 995
        syn_IN1dend_IN2.Fac = 3
        syn_IN1dend_IN2.u0 = 0.5
        nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[85](1)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[85])
        nc_IN1dend_IN2.weight[0] = 10 #5
        nc_IN1dend_IN2.delay = 0.1
        nc_IN1dend_IN2.threshold = -34
    
        syn_IN2dend_IN1 = h.DetGABAAB(subnetwork.IN_1.model.dend[25](1))
        syn_IN2dend_IN1.tau_r_GABAA = 0.71
        syn_IN2dend_IN1.tau_d_GABAA = 4.18
        syn_IN2dend_IN1.e_GABAA = -80
        syn_IN2dend_IN1.Use = 0.25
        syn_IN2dend_IN1.Dep = 995
        syn_IN2dend_IN1.Fac = 3
        syn_IN2dend_IN1.u0 = 0.5
        nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[85](1)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[85])
        nc_IN2dend_IN1.weight[0] = 10 #5
        nc_IN2dend_IN1.delay = 0.1
        nc_IN2dend_IN1.threshold = -34
        
        # add Br inputs to other INs (diffIN network arrangement)
        synBr1toIN2 = h.DetAMPANMDA(subnetwork2.IN_1.model.dend[25](1))
        synBr1toIN2.tau_r_AMPA = 0.37
        synBr1toIN2.tau_d_AMPA = 1.65
        synBr1toIN2.e = 10
        synBr1toIN2.Use = 0.48
        synBr1toIN2.Dep = 690
        synBr1toIN2.Fac = 57
        synBr1toIN2.u0 = 0.5
        nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        nc_Br1toIN2.delay = 1
        nc_Br1toIN2.weight[0] = 8 #4
        nc2_Br1toIN2 = h.NetCon(subnetwork._vecstim2, synBr1toIN2)
        nc2_Br1toIN2.delay = 1
        nc2_Br1toIN2.weight[0] = 8 #4

    
        synBr1toIN1 = h.DetAMPANMDA(subnetwork.IN_1.model.dend[85](1))
        synBr1toIN1.tau_r_AMPA = 0.36
        synBr1toIN1.tau_d_AMPA = 1.77
        synBr1toIN1.e = 10
        synBr1toIN1.Use = 0.48
        synBr1toIN1.Dep = 690
        synBr1toIN1.Fac = 57
        synBr1toIN1.u0 = 0.5
        nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        nc_Br1toIN1.delay = 1
        nc_Br1toIN1.weight[0] = 5 #2.5
        nc2_Br1toIN1 = h.NetCon(subnetwork._vecstim2, synBr1toIN1)
        nc2_Br1toIN1.delay = 1
        nc2_Br1toIN1.weight[0] = 5 #2.5
    
        synBr2toIN1 = h.DetAMPANMDA(subnetwork.IN_1.model.dend[25](1))
        synBr2toIN1.tau_r_AMPA = 0.37
        synBr2toIN1.tau_d_AMPA = 1.65
        synBr2toIN1.e = 10
        synBr2toIN1.Use = 0.48
        synBr2toIN1.Dep = 690
        synBr2toIN1.Fac = 57
        synBr2toIN1.u0 = 0.5
        nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        nc_Br2toIN1.delay = 1 
        nc_Br2toIN1.weight[0] = 8 #4
        nc2_Br2toIN1 = h.NetCon(subnetwork2._vecstim2, synBr2toIN1)
        nc2_Br2toIN1.delay = 1 
        nc2_Br2toIN1.weight[0] = 8 #4
    
        synBr2toIN2 = h.DetAMPANMDA(subnetwork2.IN_1.model.dend[85](1))
        synBr2toIN2.tau_r_AMPA = 0.36
        synBr2toIN2.tau_d_AMPA = 1.77
        synBr2toIN2.e = 10
        synBr2toIN2.Use = 0.48
        synBr2toIN2.Dep = 690
        synBr2toIN2.Fac = 57
        synBr2toIN2.u0 = 0.5
        nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        nc_Br2toIN2.delay = 1 
        nc_Br2toIN2.weight[0] = 5 #2.5
        nc2_Br2toIN2 = h.NetCon(subnetwork2._vecstim2, synBr2toIN2)
        nc2_Br2toIN2.delay = 1 
        nc2_Br2toIN2.weight[0] = 5 #2.5
        
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
        folder = '../spike_output/' + parent_dir + '/diffrun_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(x)+pickle_ext, diffIN_spike_events, allow_pickle = True)

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
    # del conversion
    # del converted_time
    # del datetimestr
    # del diffIN_spike_events
    # del diffIN_time_series
    # del diffvertstack2
    # del end
    # del endtime
    # del filename_spikes
    # del folder
    # del nc_Br1toIN1
    # del nc_Br1toIN2
    # del nc_Br2toIN1
    # del nc_Br2toIN2
    # del nc_IN1axon_IN2
    # del nc_IN1dend_IN2
    # del nc_IN2axon_IN1
    # del nc_IN2dend_IN1
    # del new_subnetwork1
    # del new_subnetwork2
    # del path_diff
    # del pickle_ext
    # del seconds_input
    # del start
    # del starttime
    # del synBr1toIN1
    # del synBr1toIN2
    # del synBr2toIN1
    # del synBr2toIN2
    # del syn_IN1axon_IN2
    # del syn_IN1dend_IN2
    # del syn_IN2axon_IN1
    # del syn_IN2dend_IN1
    # del t
    # del x
    return z
    
def samerun(z):     
    parent_dir = 'run_200711_sameopp_stp'
    
    N=4
    subnetwork = Subnetwork(0,N,100)
    subnetwork2 = Subnetwork(9,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)    
    
    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new_subnetwork2 = np.load('../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork2 = np.load('../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_same = '../spike_output/' + parent_dir + '/samerun_' + str(z) + '/'
    if not os.path.exists(path_same):
        os.mkdir(path_same)
        print("Directory " , path_same ,  " Created ")
            
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    sameIN_time_series = []
    sameIN_spike_events = []
    for y in range(0,16): #loop through different spot diameter spike trains
        # change spike inputs to network
        subnetwork.train_vec = h.Vector(new_subnetwork1[y,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[y,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        subnetwork2.train_vec = h.Vector(new_subnetwork2[y,])
        subnetwork2._vecstim.play(subnetwork2.train_vec)
        subnetwork2.train_vec2 = h.Vector(new2_subnetwork2[y,])
        subnetwork2._vecstim2.play(subnetwork2.train_vec2)
        # here is where I connect the two subnetworks with each other
        # connect two subnetworks by adding IN-IN axons and dendrites
        syn_IN1axon_IN2 = h.DetGABAAB(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1axon_IN2.tau_r_GABAA = 0.71
        syn_IN1axon_IN2.tau_d_GABAA = 4.18
        syn_IN1axon_IN2.e_GABAA = -80
        syn_IN1axon_IN2.Use = 0.25
        syn_IN1axon_IN2.Dep = 995
        syn_IN1axon_IN2.Fac = 3
        syn_IN1axon_IN2.u0 = 0.5
        nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
        nc_IN1axon_IN2.weight[0] = 10 #5
        nc_IN1axon_IN2.delay = 1
    
        syn_IN2axon_IN1 = h.DetGABAAB(subnetwork.IN_1.model.dend[25](1))
        syn_IN2axon_IN1.tau_r_GABAA = 0.71
        syn_IN2axon_IN1.tau_d_GABAA = 4.18
        syn_IN2axon_IN1.e_GABAA = -80
        syn_IN2axon_IN1.Use = 0.25
        syn_IN2axon_IN1.Dep = 995
        syn_IN2axon_IN1.Fac = 3
        syn_IN2axon_IN1.u0 = 0.5
        nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
        nc_IN2axon_IN1.weight[0] = 10 #5
        nc_IN2axon_IN1.delay = 1
    
        syn_IN1dend_IN2 = h.DetGABAAB(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1dend_IN2.tau_r_GABAA = 0.71
        syn_IN1dend_IN2.tau_d_GABAA = 4.18
        syn_IN1dend_IN2.e_GABAA = -80
        syn_IN1dend_IN2.Use = 0.25
        syn_IN1dend_IN2.Dep = 995
        syn_IN1dend_IN2.Fac = 3
        syn_IN1dend_IN2.u0 = 0.5
        nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[85](1)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[85])
        nc_IN1dend_IN2.weight[0] = 10 #5
        nc_IN1dend_IN2.delay = 0.1
        nc_IN1dend_IN2.threshold = -34
    
        syn_IN2dend_IN1 = h.DetGABAAB(subnetwork.IN_1.model.dend[25](1))
        syn_IN2dend_IN1.tau_r_GABAA = 0.71
        syn_IN2dend_IN1.tau_d_GABAA = 4.18
        syn_IN2dend_IN1.e_GABAA = -80
        syn_IN2dend_IN1.Use = 0.25
        syn_IN2dend_IN1.Dep = 995
        syn_IN2dend_IN1.Fac = 3
        syn_IN2dend_IN1.u0 = 0.5
        nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[85](1)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[85])
        nc_IN2dend_IN1.weight[0] = 10 #5
        nc_IN2dend_IN1.delay = 0.1
        nc_IN2dend_IN1.threshold = -34
       
        # add Br inputs to own IN (sameIN network arrangement)
        synBr1toIN1 = h.DetAMPANMDA(subnetwork.IN_1.model.dend[25](1))
        synBr1toIN1.tau_r_AMPA = 0.37
        synBr1toIN1.tau_d_AMPA = 1.65
        synBr1toIN1.e = 10 
        synBr1toIN1.Use = 0.48
        synBr1toIN1.Dep = 690
        synBr1toIN1.Fac = 57
        synBr1toIN1.u0 = 0.5
        nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        nc_Br1toIN1.delay = 1 
        nc_Br1toIN1.weight[0] = 8 #4
        nc2_Br1toIN1 = h.NetCon(subnetwork._vecstim2, synBr1toIN1)
        nc2_Br1toIN1.delay = 1 
        nc2_Br1toIN1.weight[0] = 8 #4
    
        synBr1toIN2 = h.DetAMPANMDA(subnetwork2.IN_1.model.dend[85](1))
        synBr1toIN2.tau_r_AMPA = 0.36
        synBr1toIN2.tau_d_AMPA = 1.77
        synBr1toIN2.e = 10
        synBr1toIN2.Use = 0.48
        synBr1toIN2.Dep = 690
        synBr1toIN2.Fac = 57
        synBr1toIN2.u0 = 0.5
        nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        nc_Br1toIN2.delay = 1 
        nc_Br1toIN2.weight[0] = 5 #2.5
        nc2_Br1toIN2 = h.NetCon(subnetwork._vecstim2, synBr1toIN2)
        nc2_Br1toIN2.delay = 1 
        nc2_Br1toIN2.weight[0] = 5 #2.5
    
        synBr2toIN2 = h.DetAMPANMDA(subnetwork2.IN_1.model.dend[25](1))
        synBr2toIN2.tau_r_AMPA = 0.37
        synBr2toIN2.tau_d_AMPA = 1.65
        synBr2toIN2.e = 10
        synBr2toIN2.Use = 0.48
        synBr2toIN2.Dep = 690
        synBr2toIN2.Fac = 57
        synBr2toIN2.u0 = 0.5
        nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        nc_Br2toIN2.delay = 1 
        nc_Br2toIN2.weight[0] = 8 #4
        nc2_Br2toIN2 = h.NetCon(subnetwork2._vecstim2, synBr2toIN2)
        nc2_Br2toIN2.delay = 1 
        nc2_Br2toIN2.weight[0] = 8 #4
    
        synBr2toIN1 = h.DetAMPANMDA(subnetwork.IN_1.model.dend[85](1))
        synBr2toIN1.tau_r_AMPA = 0.36
        synBr2toIN1.tau_d_AMPA = 1.77
        synBr2toIN1.e = 10
        synBr2toIN1.Use = 0.48
        synBr2toIN1.Dep = 690
        synBr2toIN1.Fac = 57
        synBr2toIN1.u0 = 0.5
        nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        nc_Br2toIN1.delay = 1 
        nc_Br2toIN1.weight[0] = 5 #2.5
        nc2_Br2toIN1 = h.NetCon(subnetwork2._vecstim2, synBr2toIN1)
        nc2_Br2toIN1.delay = 1 
        nc2_Br2toIN1.weight[0] = 5 #2.5
        
        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()
        #<3  samevertstack = []
        #<3  samevertstack = np.vstack((
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
        samevertstack2 = []
        samevertstack2 = [
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
        
        #<3 sameIN_time_series.append(samevertstack)
        #<3 sameIN_spike_events.append(samevertstack2)
        #<3 sameIN_spike_events = np.array(sameIN_spike_events)
        sameIN_spike_events = np.array(samevertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        #<3 filename_timeseries = 'sameIN_time_series'
        filename_spikes = 'sameIN_spike_events' 
        #<3 txt_ext = '.txt'
        pickle_ext = '.npy'
        folder = '../spike_output/' + parent_dir + '/samerun_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(y)+pickle_ext, sameIN_spike_events, allow_pickle = True)
            
    #<3 sameIN_time_series = np.array(sameIN_time_series)
    
    #<3 print(sameIN_time_series.shape)
    #<3 print(sameIN_spike_events.shape)
    
    # Write the arrays to disk
    
    #<3 with open(folder+filename_timeseries+datetimestr+txt_ext, 'w') as outfile:
    #<3     outfile.write('# Array shape: {0}\n'.format(sameIN_time_series.shape))
    #<3     for data_slice in sameIN_time_series:
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

    # del conversion
    # del converted_time
    # del datetimestr
    # del end
    # del endtime
    # del filename_spikes
    # del folder
    # del nc_Br1toIN1
    # del nc_Br1toIN2
    # del nc_Br2toIN1
    # del nc_Br2toIN2
    # del nc_IN1axon_IN2
    # del nc_IN1dend_IN2
    # del nc_IN2axon_IN1
    # del nc_IN2dend_IN1
    # del new_subnetwork1
    # del new_subnetwork2
    # del path_same
    # del pickle_ext
    # del sameIN_spike_events
    # del sameIN_time_series
    # del samevertstack2
    # del seconds_input
    # del start
    # del starttime
    # del synBr1toIN1
    # del synBr1toIN2
    # del synBr2toIN1
    # del synBr2toIN2
    # del syn_IN1axon_IN2
    # del syn_IN1dend_IN2
    # del syn_IN2axon_IN1
    # del syn_IN2dend_IN1
    # del t
    # del y 
    return z

def makeInputs(z):
    parent_dir = 'run_200711_sameopp_stp/'


    timenow = datetime.datetime.now().time()
    print('now reading: '+str(z) + '  ',timenow)

    poisson_whisker.createSubnetworkInputs(z,parent_dir)

    return z
    
if __name__ == '__main__':  
    pool = multiprocessing.Pool(20)
    z = [1,2,3,4,5,6,7,8,9,10]
    # result1 = pool.map(makeInputs,z)
    # parent_dir = 'run_200711_sameopp_stp/'
    # convert_inputs.convertinputs(parent_dir)
    result2 = pool.map_async(diffrun,z)
    result3 = pool.map_async(samerun,z)
    print(result2.get())
    print(result3.get())
    freq = 440
    dur = 1000
    winsound.Beep(freq,dur)
