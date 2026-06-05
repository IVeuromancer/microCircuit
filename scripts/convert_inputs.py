import glob
import numpy as np

def convertinputs(parent_dir_str):

    spike_times_array_sub1 = []
    spike_times_array_sub2 = []
    spike_times = []
    to_be_averaged = []

    path = '../spike_input/' +parent_dir_str+'/Br1/subnetwork1_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub1 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub1_new = []
    spike_times_array_sub1_new = spike_times_array_sub1[:,16:] 

    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'Br1/edited/subnetwork1_input' + '_' + trial_str + '.npy', spike_times_array_sub1_new[idxnum], allow_pickle = True)

    path = '../spike_input/' +parent_dir_str +'/Br1/subnetwork2_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub2 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub2_new = []
    spike_times_array_sub2_new = spike_times_array_sub2[:,16:32]
    spike_times_array_sub2_new_new = np.hstack((spike_times_array_sub2_new, spike_times_array_sub2[:,0:16]))

    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'Br1/edited/subnetwork2_input' + '_' + trial_str + '.npy', spike_times_array_sub2_new_new[idxnum], allow_pickle = True)

    path = '../spike_input/' +parent_dir_str+'/Br2/subnetwork1_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub1 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub1_new = []
    spike_times_array_sub1_new = spike_times_array_sub1[:,16:]

    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'Br2/edited/subnetwork1_input' + '_' + trial_str + '.npy', spike_times_array_sub1_new[idxnum], allow_pickle = True)

    path = '../spike_input/' +parent_dir_str+'/Br2/subnetwork2_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub2 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub2_new = []
    spike_times_array_sub2_new = spike_times_array_sub2[:,16:32]
    spike_times_array_sub2_new_new = np.hstack((spike_times_array_sub2_new, spike_times_array_sub2[:,0:16]))


    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'Br2/edited/subnetwork2_input' + '_' + trial_str + '.npy', spike_times_array_sub2_new_new[idxnum], allow_pickle = True)


    path = '../spike_input/' +parent_dir_str+'/extra_input1/subnetworkall_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub1 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub1_new = []
    spike_times_array_sub1_new = spike_times_array_sub1[:,16:]

    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'extra_input1/edited/subnetworkall_input' + '_' + trial_str + '.npy', spike_times_array_sub1_new[idxnum], allow_pickle = True)


    path = '../spike_input/' +parent_dir_str+'/extra_input2/subnetworkall_input*'
    for filename in glob.glob(path):
        spikes = np.load(filename, allow_pickle = True)
        spike_times.append(spikes) 
    spike_times_array_sub1 = np.array(spike_times)
    spike_times = []
    spike_times_array_sub1_new = []
    spike_times_array_sub1_new = spike_times_array_sub1[:,16:]

    for x in range(1,11):
        trial_str = str(x)
        idxnum = x-1
        np.save('../spike_input/' + parent_dir_str+'/'+'extra_input2/edited/subnetworkall_input' + '_' + trial_str + '.npy', spike_times_array_sub1_new[idxnum], allow_pickle = True)