import poisson
import os
import diffrun_200525
import samerun_200526
import winsound

parent_dir = 'run_200525'

for z in range(5,6):
    poisson.createSubnetworkInputs(trial = z) 
    
    path_diff = '../../spike_output/' + parent_dir + '/diffrun_' + str(z)
    if not os.path.exists(path_diff):
        os.mkdir(path_diff)
        print("Directory " , path_diff ,  " Created ")
 
    diffrun_200525.runDiff(parent_dir, z)
    
    path_same = '../../spike_output/' + parent_dir + '/samerun_' + str(z)
    if not os.path.exists(path_same):
        os.mkdir(path_same)
        print("Directory " , path_same ,  " Created ")
    
    samerun_200526.runSame(parent_dir, z)

freq = 440
dur = 3000
winsound.Beep(freq,dur)