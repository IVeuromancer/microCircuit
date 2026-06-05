from neuron import h
from subnetwork_200620_heiberg import Subnetwork
from subnetwork_200620_heiberg_no_triadic import Subnetwork_no_triadic
from subnetwork_200620_heiberg_no_axonal import Subnetwork_no_axonal
from subnetwork_200620_heiberg_no_inhib import Subnetwork_no_inhib
import poisson
import numpy as np
import time
import datetime
import winsound
import os
import multiprocessing
from multiprocessing import freeze_support, Pool

def runHeiberg(z):
    parent_dir = 'run_200704_heiberg'

    N=4
    subnetwork = Subnetwork(0,N,100)

    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): 

        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()

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
        np.array(subnetwork.TCcells[3].spike_times)]
        
        diffIN_spike_events = np.array(diffvertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        filename_spikes = 'spike_events' 
        pickle_ext = '.npy'
        folder = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(x)+pickle_ext, diffIN_spike_events, allow_pickle = True)
    
    endtime = datetime.datetime.now().time()
    print ("end time :", endtime)
    end = time.time()
    seconds_input = end - start
    conversion = datetime.timedelta(seconds=seconds_input)
    converted_time = str(conversion)
    # print(seconds_input)
    print(converted_time)
    return z

def runHeiberg_notriadic(z):
    parent_dir = 'run_200704_heiberg_no_triadic'

    N=4
    subnetwork = Subnetwork_no_triadic(0,N,100)

    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): 

        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()

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
        np.array(subnetwork.TCcells[3].spike_times)]
        
        diffIN_spike_events = np.array(diffvertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        filename_spikes = 'spike_events' 
        pickle_ext = '.npy'
        folder = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(x)+pickle_ext, diffIN_spike_events, allow_pickle = True)
    
    endtime = datetime.datetime.now().time()
    print ("end time :", endtime)
    end = time.time()
    seconds_input = end - start
    conversion = datetime.timedelta(seconds=seconds_input)
    converted_time = str(conversion)
    # print(seconds_input)
    print(converted_time)
    return z

def runHeiberg_noaxonal(z):
    parent_dir = 'run_200704_heiberg_no_axonal'

    N=4
    subnetwork = Subnetwork_no_axonal(0,N,100)

    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): 

        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()

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
        np.array(subnetwork.TCcells[3].spike_times)]
        
        diffIN_spike_events = np.array(diffvertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        filename_spikes = 'spike_events' 
        pickle_ext = '.npy'
        folder = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(x)+pickle_ext, diffIN_spike_events, allow_pickle = True)
    
    endtime = datetime.datetime.now().time()
    print ("end time :", endtime)
    end = time.time()
    seconds_input = end - start
    conversion = datetime.timedelta(seconds=seconds_input)
    converted_time = str(conversion)
    # print(seconds_input)
    print(converted_time)
    return z

def runHeiberg_noinhib(z):
    parent_dir = 'run_200704_heiberg_no_inhib'

    N=4
    subnetwork = Subnetwork_no_inhib(0,N,100)

    new_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br1/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)
    new2_subnetwork1 = np.load('../spike_input/'+parent_dir+'/'+'Br2/subnetwork1_input' + '_' + str(z) + '.npy', allow_pickle = True)

    path_diff = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
    
    start = time.time()
    starttime = datetime.datetime.now().time()
    print ("start time :", starttime)
 
    diffIN_time_series = []
    diffIN_spike_events = []
    for x in range(0,len(new_subnetwork1)): 

        subnetwork.train_vec = h.Vector(new_subnetwork1[x,])
        subnetwork._vecstim.play(subnetwork.train_vec)
        subnetwork.train_vec2 = h.Vector(new2_subnetwork1[x,])
        subnetwork._vecstim2.play(subnetwork.train_vec2)

        t = h.Vector().record(h._ref_t)
        h.v_init = -78
        h.celsius = 34
        h.tstop = 1000
        h.init()
        h.run()

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
        np.array(subnetwork.TCcells[3].spike_times)]
        
        diffIN_spike_events = np.array(diffvertstack2)
        datetimestr = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
        filename_spikes = 'spike_events' 
        pickle_ext = '.npy'
        folder = '../spike_output/' + parent_dir + '/run_' + str(z) + '/'
        np.save(folder+filename_spikes+datetimestr+'-'+str(x)+pickle_ext, diffIN_spike_events, allow_pickle = True)
    
    endtime = datetime.datetime.now().time()
    print ("end time :", endtime)
    end = time.time()
    seconds_input = end - start
    conversion = datetime.timedelta(seconds=seconds_input)
    converted_time = str(conversion)
    # print(seconds_input)
    print(converted_time)
    return z

def makeInputs(z):
    parent_dir = 'run_200627_heiberg/'


    timenow = datetime.datetime.now().time()
    print('now reading: '+str(z) + '  ',timenow)

    poisson.createSubnetworkInputs(z,parent_dir)

    return z



    
if __name__ == '__main__':  
    freeze_support()

    z = [1,2,3,4,5,6,7,8,9,10]
    pool = multiprocessing.Pool(40)
    # result1 = pool.map(makeInputs,z)

    result2 = pool.map_async(runHeiberg,z)
    result3 = pool.map_async(runHeiberg_noaxonal,z)
    result4 = pool.map_async(runHeiberg_notriadic,z)
    result5 = pool.map_async(runHeiberg_noinhib,z)
    print(result2.get())
    print(result3.get())
    print(result4.get())
    print(result5.get())

    freq = 440
    dur = 1000
    winsound.Beep(freq,dur)
    
    














