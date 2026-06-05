from neuron import h
from subnetwork_200619_no_inhib import Subnetwork_no_inhib
import numpy as np
import time
import datetime
import winsound
import os
from distutils.dir_util import copy_tree

def copyinputs(old_parent_dir):
    datetimestr = datetime.datetime.now().strftime("_%y%m%d")  
    parent_dir = 'run' + datetimestr + '_sameopp_no_inhib'
    old_dir = '../../spike_input/' + old_parent_dir
    new_dir = '../../spike_input/' + parent_dir  
    copy_tree(old_dir, new_dir)

def diffrun(z):    

    #datetimestr = datetime.datetime.now().strftime("_%y%m%d")  
    datetimestr = '_200914'
    parent_dir = 'run' +datetimestr + '_sameopp_no_inhib' 

    # path_diff = '../../spike_output/' + parent_dir 
    # if not os.path.exists(path_diff):
    #     os.mkdir(path_diff)
    #     print("Directory " , path_diff ,  " Created ")

    N=4
    subnetwork = Subnetwork_no_inhib(0,N,100)
    subnetwork2 = Subnetwork_no_inhib(9,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)

    new_subnetwork1 = np.load('../../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new_subnetwork2 = np.load('../../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork2 = np.load('../../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new3_subnetworkall = np.load('../../spike_input/'+parent_dir+'/'+'extra_input1/edited/subnetworkall_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new4_subnetworkall = np.load('../../spike_input/'+parent_dir+'/'+'extra_input2/edited/subnetworkall_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../../spike_output/' + parent_dir + '/diffrun_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): #loop through different spot diameter spike trains
        # change spike inputs to network
        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork2.train_vec = h.Vector(new_subnetwork2[x,])
        subnetwork2._vecstim.play(subnetwork2.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)
        subnetwork2.train_vec2 = h.Vector(new2_subnetwork2[x,])
        subnetwork2._vecstim2.play(subnetwork2.train_vec2)

        subnetwork.train_vec3 = h.Vector(new3_subnetworkall[x,])
        subnetwork._vecstim3.play(subnetwork.train_vec3)
        subnetwork2.train_vec3 = h.Vector(new3_subnetworkall[x,])
        subnetwork2._vecstim3.play(subnetwork2.train_vec3)
        subnetwork.train_vec4 = h.Vector(new4_subnetworkall[x,])
        subnetwork._vecstim4.play(subnetwork.train_vec4)
        subnetwork2.train_vec4 = h.Vector(new4_subnetworkall[x,])
        subnetwork2._vecstim4.play(subnetwork2.train_vec4)


        # here is where I connect the two subnetworks with each other
        # connect two subnetworks by adding IN-IN axons and dendrites
        syn_IN1axon_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1axon_IN2.tau1 = 0.71
        syn_IN1axon_IN2.tau2 = 4.18
        syn_IN1axon_IN2.e = -80
        nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
        nc_IN1axon_IN2.weight[0] = 0.005
        nc_IN1axon_IN2.delay = 1
    
        syn_IN2axon_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2axon_IN1.tau1 = 0.71
        syn_IN2axon_IN1.tau2 = 4.18
        syn_IN2axon_IN1.e = -80
        nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
        nc_IN2axon_IN1.weight[0] = 0.005
        nc_IN2axon_IN1.delay = 1
    
        syn_IN1dend_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1dend_IN2.tau1 = 0.71
        syn_IN1dend_IN2.tau2 = 4.18
        syn_IN1dend_IN2.e = -80
        nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[85](1)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[85])
        nc_IN1dend_IN2.weight[0] = 0.005
        nc_IN1dend_IN2.delay = 0.1
        nc_IN1dend_IN2.threshold = -34
    
        syn_IN2dend_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2dend_IN1.tau1 = 0.71
        syn_IN2dend_IN1.tau2 = 4.18
        syn_IN2dend_IN1.e = -80
        nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[85](1)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[85])
        nc_IN2dend_IN1.weight[0] = 0.005
        nc_IN2dend_IN1.delay = 0.1
        nc_IN2dend_IN1.threshold = -34
        
        # add Br inputs to other INs (diffIN network arrangement)

        #temp! extra input
        syn_extra = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.5))
        syn_extra.tau1 = 0.37
        syn_extra.tau2 = 1.65
        syn_extra.e = 10
        nc_extra = h.NetCon(subnetwork._vecstim3, syn_extra)
        nc_extra.delay = 1
        nc_extra.weight[0] = 0.004
        nc2_extra = h.NetCon(subnetwork._vecstim4, syn_extra)
        nc2_extra.delay = 1
        nc2_extra.weight[0] = 0.004

        syn2_extra = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.5))
        syn2_extra.tau1 = 0.37
        syn2_extra.tau2 = 1.65
        syn2_extra.e = 10
        nc2_extra = h.NetCon(subnetwork2._vecstim3, syn2_extra)
        nc2_extra.delay = 1
        nc2_extra.weight[0] = 0.004
        nc22_extra = h.NetCon(subnetwork2._vecstim4, syn2_extra)
        nc22_extra.delay = 1
        nc22_extra.weight[0] = 0.004
        #above is temp, remove later!


        synBr1toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        synBr1toIN2.tau1 = 0.37
        synBr1toIN2.tau2 = 1.65
        synBr1toIN2.e = 10
        nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        nc_Br1toIN2.delay = 1
        nc_Br1toIN2.weight[0] = 0.004
        nc2_Br1toIN2 = h.NetCon(subnetwork._vecstim2, synBr1toIN2)
        nc2_Br1toIN2.delay = 1
        nc2_Br1toIN2.weight[0] = 0.004
        #extra inputs
        # ncsecond_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        # ncsecond_Br1toIN2.delay = 1
        # ncsecond_Br1toIN2.weight[0] = 0.004
        # ncsecond2_Br1toIN2 = h.NetCon(subnetwork._vecstim2, synBr1toIN2)
        # ncsecond2_Br1toIN2.delay = 1
        # ncsecond2_Br1toIN2.weight[0] = 0.004

    
        synBr1toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[85](1))
        synBr1toIN1.tau1 = 0.36
        synBr1toIN1.tau2 = 1.77
        synBr1toIN1.e = 10
        nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        nc_Br1toIN1.delay = 1
        nc_Br1toIN1.weight[0] = 0.0025
        nc2_Br1toIN1 = h.NetCon(subnetwork._vecstim2, synBr1toIN1)
        nc2_Br1toIN1.delay = 1
        nc2_Br1toIN1.weight[0] = 0.0025
    
        synBr2toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        synBr2toIN1.tau1 = 0.37
        synBr2toIN1.tau2 = 1.65
        synBr2toIN1.e = 10
        nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        nc_Br2toIN1.delay = 1 
        nc_Br2toIN1.weight[0] = 0.004
        nc2_Br2toIN1 = h.NetCon(subnetwork2._vecstim2, synBr2toIN1)
        nc2_Br2toIN1.delay = 1 
        nc2_Br2toIN1.weight[0] = 0.004
        #extra inputs
        ncsecond_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        ncsecond_Br2toIN1.delay = 1 
        ncsecond_Br2toIN1.weight[0] = 0.004
        ncsecond2_Br2toIN1 = h.NetCon(subnetwork2._vecstim2, synBr2toIN1)
        ncsecond2_Br2toIN1.delay = 1 
        ncsecond2_Br2toIN1.weight[0] = 0.004
    
        synBr2toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[85](1))
        synBr2toIN2.tau1 = 0.36
        synBr2toIN2.tau2 = 1.77
        synBr2toIN2.e = 10
        nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        nc_Br2toIN2.delay = 1 
        nc_Br2toIN2.weight[0] = 0.0025
        nc2_Br2toIN2 = h.NetCon(subnetwork2._vecstim2, synBr2toIN2)
        nc2_Br2toIN2.delay = 1 
        nc2_Br2toIN2.weight[0] = 0.0025
        
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
        folder = '../../spike_output/' + parent_dir + '/diffrun_' + str(z) + '/'
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
    # del ncsecond_Br1toIN2
    # del new_subnetwork1
    # del new_subnetwork2
    # del new2_subnetwork1
    # del new2_subnetwork2
    # del path_diff
    # del pickle_ext
    # del seconds_input
    # del start
    # del starttime
    # del subnetwork
    # del subnetwork2
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
    #datetimestr = datetime.datetime.now().strftime("_%y%m%d")  
    datetimestr = '_200914'
    parent_dir = 'run' +datetimestr + '_sameopp_no_inhib' 

    # path_same = '../../spike_output/' + parent_dir 
    # if not os.path.exists(path_same):
    #     os.mkdir(path_same)
    #     print("Directory " , path_same ,  " Created ")
    
    N=4
    subnetwork = Subnetwork_no_inhib(0,N,100)
    subnetwork2 = Subnetwork_no_inhib(9,N,100) #gidStart should be number of gids (3 for IN + 2 for TC)    
    
    new_subnetwork1 = np.load('../../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new_subnetwork2 = np.load('../../spike_input/'+parent_dir+'/'+'Br1/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork2 = np.load('../../spike_input/'+parent_dir+'/'+'Br2/edited/subnetwork2_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new3_subnetworkall = np.load('../../spike_input/'+parent_dir+'/'+'extra_input1/edited/subnetworkall_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new4_subnetworkall = np.load('../../spike_input/'+parent_dir+'/'+'extra_input2/edited/subnetworkall_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_same = '../../spike_output/' + parent_dir + '/samerun_' + str(z) + '/'
    if not os.path.exists(path_same):
        os.mkdir(path_same)
        print("Directory " , path_same ,  " Created ")
            
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    sameIN_time_series = []
    sameIN_spike_events = []
    for y in range(0,len(new_subnetwork1)): #loop through different spot diameter spike trains
        # change spike inputs to network
        subnetwork.train_vec = h.Vector(new_subnetwork1[y,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork2.train_vec = h.Vector(new_subnetwork2[y,])
        subnetwork2._vecstim.play(subnetwork2.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[y,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)
        subnetwork2.train_vec2 = h.Vector(new2_subnetwork2[y,])
        subnetwork2._vecstim2.play(subnetwork2.train_vec2)

        subnetwork.train_vec3 = h.Vector(new3_subnetworkall[y,])
        subnetwork._vecstim3.play(subnetwork.train_vec3)
        subnetwork2.train_vec3 = h.Vector(new3_subnetworkall[y,])
        subnetwork2._vecstim3.play(subnetwork2.train_vec3)
        subnetwork.train_vec4 = h.Vector(new4_subnetworkall[y,])
        subnetwork._vecstim4.play(subnetwork.train_vec4)
        subnetwork2.train_vec4 = h.Vector(new4_subnetworkall[y,])
        subnetwork2._vecstim4.play(subnetwork2.train_vec4)
        # here is where I connect the two subnetworks with each other
        # connect two subnetworks by adding IN-IN axons and dendrites
        syn_IN1axon_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1axon_IN2.tau1 = 0.71
        syn_IN1axon_IN2.tau2 = 4.18
        syn_IN1axon_IN2.e = -80
        nc_IN1axon_IN2 = h.NetCon(subnetwork.IN_1.model.soma[0](0.5)._ref_v, syn_IN1axon_IN2, sec=subnetwork.IN_1.model.soma[0])
        nc_IN1axon_IN2.weight[0] = 0.005
        nc_IN1axon_IN2.delay = 1
    
        syn_IN2axon_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2axon_IN1.tau1 = 0.71
        syn_IN2axon_IN1.tau2 = 4.18
        syn_IN2axon_IN1.e = -80
        nc_IN2axon_IN1 = h.NetCon(subnetwork2.IN_1.model.soma[0](0.5)._ref_v, syn_IN2axon_IN1, sec=subnetwork2.IN_1.model.soma[0])
        nc_IN2axon_IN1.weight[0] = 0.005
        nc_IN2axon_IN1.delay = 1
    
        syn_IN1dend_IN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        syn_IN1dend_IN2.tau1 = 0.71
        syn_IN1dend_IN2.tau2 = 4.18
        syn_IN1dend_IN2.e = -80
        nc_IN1dend_IN2 = h.NetCon(subnetwork.IN_1.model.dend[85](1)._ref_v, syn_IN1dend_IN2, sec=subnetwork.IN_1.model.dend[85])
        nc_IN1dend_IN2.weight[0] = 0.005
        nc_IN1dend_IN2.delay = 0.1
        nc_IN1dend_IN2.threshold = -34
    
        syn_IN2dend_IN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        syn_IN2dend_IN1.tau1 = 0.71
        syn_IN2dend_IN1.tau2 = 4.18
        syn_IN2dend_IN1.e = -80
        nc_IN2dend_IN1 = h.NetCon(subnetwork2.IN_1.model.dend[85](1)._ref_v, syn_IN2dend_IN1, sec=subnetwork2.IN_1.model.dend[85])
        nc_IN2dend_IN1.weight[0] = 0.005
        nc_IN2dend_IN1.delay = 0.1
        nc_IN2dend_IN1.threshold = -34
       
        # add Br inputs to own IN (sameIN network arrangement)

        #temp! extra input
        syn_extra = h.Exp2Syn(subnetwork.IN_1.model.dend[0](0.5))
        syn_extra.tau1 = 0.37
        syn_extra.tau2 = 1.65
        syn_extra.e = 10
        nc_extra = h.NetCon(subnetwork._vecstim3, syn_extra)
        nc_extra.delay = 1
        nc_extra.weight[0] = 0.004
        nc2_extra = h.NetCon(subnetwork._vecstim4, syn_extra)
        nc2_extra.delay = 1
        nc2_extra.weight[0] = 0.004

        syn2_extra = h.Exp2Syn(subnetwork2.IN_1.model.dend[0](0.5))
        syn2_extra.tau1 = 0.37
        syn2_extra.tau2 = 1.65
        syn2_extra.e = 10
        nc2_extra = h.NetCon(subnetwork2._vecstim3, syn2_extra)
        nc2_extra.delay = 1
        nc2_extra.weight[0] = 0.004
        nc22_extra = h.NetCon(subnetwork2._vecstim4, syn2_extra)
        nc22_extra.delay = 1
        nc22_extra.weight[0] = 0.004
        #above is temp, remove later!

        synBr1toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[25](1))
        synBr1toIN1.tau1 = 0.37
        synBr1toIN1.tau2 = 1.65
        synBr1toIN1.e = 10 
        nc_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        nc_Br1toIN1.delay = 1 
        nc_Br1toIN1.weight[0] = 0.004
        nc2_Br1toIN1 = h.NetCon(subnetwork._vecstim2, synBr1toIN1)
        nc2_Br1toIN1.delay = 1 
        nc2_Br1toIN1.weight[0] = 0.004
        #extra input
        # ncsecond_Br1toIN1 = h.NetCon(subnetwork._vecstim, synBr1toIN1)
        # ncsecond_Br1toIN1.delay = 1 
        # ncsecond_Br1toIN1.weight[0] = 0.004
        # ncsecond2_Br1toIN1 = h.NetCon(subnetwork._vecstim2, synBr1toIN1)
        # ncsecond2_Br1toIN1.delay = 1 
        # ncsecond2_Br1toIN1.weight[0] = 0.004
    
        synBr1toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[85](1))
        synBr1toIN2.tau1 = 0.36
        synBr1toIN2.tau2 = 1.77
        synBr1toIN2.e = 10
        nc_Br1toIN2 = h.NetCon(subnetwork._vecstim, synBr1toIN2)
        nc_Br1toIN2.delay = 1 
        nc_Br1toIN2.weight[0] = 0.0025
        nc2_Br1toIN2 = h.NetCon(subnetwork._vecstim2, synBr1toIN2)
        nc2_Br1toIN2.delay = 1 
        nc2_Br1toIN2.weight[0] = 0.0025
    
        synBr2toIN2 = h.Exp2Syn(subnetwork2.IN_1.model.dend[25](1))
        synBr2toIN2.tau1 = 0.37
        synBr2toIN2.tau2 = 1.65
        synBr2toIN2.e = 10
        nc_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        nc_Br2toIN2.delay = 1 
        nc_Br2toIN2.weight[0] = 0.004
        nc2_Br2toIN2 = h.NetCon(subnetwork2._vecstim2, synBr2toIN2)
        nc2_Br2toIN2.delay = 1 
        nc2_Br2toIN2.weight[0] = 0.004
        #extra input
        ncsecond_Br2toIN2 = h.NetCon(subnetwork2._vecstim, synBr2toIN2)
        ncsecond_Br2toIN2.delay = 1 
        ncsecond_Br2toIN2.weight[0] = 0.004
        ncsecond2_Br2toIN2 = h.NetCon(subnetwork2._vecstim2, synBr2toIN2)
        ncsecond2_Br2toIN2.delay = 1 
        ncsecond2_Br2toIN2.weight[0] = 0.004
    
        synBr2toIN1 = h.Exp2Syn(subnetwork.IN_1.model.dend[85](1))
        synBr2toIN1.tau1 = 0.36
        synBr2toIN1.tau2 = 1.77
        synBr2toIN1.e = 10
        nc_Br2toIN1 = h.NetCon(subnetwork2._vecstim, synBr2toIN1)
        nc_Br2toIN1.delay = 1 
        nc_Br2toIN1.weight[0] = 0.0025
        nc2_Br2toIN1 = h.NetCon(subnetwork2._vecstim2, synBr2toIN1)
        nc2_Br2toIN1.delay = 1 
        nc2_Br2toIN1.weight[0] = 0.0025
        
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
        folder = '../../spike_output/' + parent_dir + '/samerun_' + str(z) + '/'
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
    # del ncsecond_Br1toIN1
    # del new_subnetwork1
    # del new_subnetwork2
    # del new2_subnetwork1
    # del new2_subnetwork2
    # del path_same
    # del pickle_ext
    # del sameIN_spike_events
    # del sameIN_time_series
    # del samevertstack2
    # del seconds_input
    # del start
    # del starttime
    # del subnetwork
    # del subnetwork2
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














