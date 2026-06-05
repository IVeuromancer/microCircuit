from neuron import h
from cell_200619 import IN
from cell_200619 import TC

class Subnetwork_no_triadic:
    """
    A network of TC and IN cells where IN makes an
    axonal and dendritic synapse onto all other TCs and all cells 
    receive Br input
    
    """
    def __init__(self, gidStart=0, N=2, r=100):
        """
        :param N: Number of cells. 2 by default
        :param r: radius of the network. 100um by default

        """
        self._N = N
        self.gidStart = gidStart #this is to offset gids every time instantiate a subnetwork (e.g. in subnetwork1, IN has gids 0,1,2 and TCs have gids 3,4. Then subnetwork2 must start with gid 5)
        
        # self.ncBrs = []
        
        # Add Br spikes of nonstationary poisson process. from poisson.py
        self._vecstim = h.VecStim()
        self._vecstim2 = h.VecStim()
        self._vecstim3 = h.VecStim()
        self._vecstim4 = h.VecStim()
        self.ncBrsPs = []
        self.nc2BrsPs = []
        
        # create TCs and IN class objects and add synapses from IN to TC
        self._create_cells(r)
        
        # create NetCon from simple spikes and poisson process (external inputs) to TC
        self._add_IN_syn()
        self._add_Br_Ps_syn()
        
        # create NetCon from Br for IN targets (prox and dist denrites)
        self.ncBr1 = h.NetCon(self._vecstim, self.IN_1.syn_Brprox)
        self.ncBr1.weight[0] = 0.004
        self.ncBr1.delay = 1
        self.nc2Br1 = h.NetCon(self._vecstim2, self.IN_1.syn_Brprox)
        self.nc2Br1.weight[0] = 0.004
        self.nc2Br1.delay = 1
        # self.ncsecondBr1 = h.NetCon(self._vecstim, self.IN_1.syn_Brprox)
        # self.ncsecondBr1.weight[0] = 0.004
        # self.ncsecondBr1.delay = 1
        # self.ncsecond2Br1 = h.NetCon(self._vecstim2, self.IN_1.syn_Brprox)
        # self.ncsecond2Br1.weight[0] = 0.004
        # self.ncsecond2Br1.delay = 1
        
        # this is from other Br input (one of the arrangements)
        # self.ncBr2 = h.NetCon(self._vecstim, self.IN_1.syn_Brprox2)
        # self.ncBr2.weight[0] =0.004
        # self.ncBr2.delay = 1
        
        self.ncBr3 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist)
        self.ncBr3.weight[0] =0.0025 #changed from 0.002
        self.ncBr3.delay = 1
        # self.nc2Br3 = h.NetCon(self._vecstim2, self.IN_1.syn_Brdist)
        # self.nc2Br3.weight[0] =0.0025 #changed from 0.002
        # self.nc2Br3.delay = 1
        
        # this is one of the arrangements
        # self.ncBr4 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist2)
        # self.ncBr4.weight[0] =0.0025
        # self.ncBr4.delay = 1
        
        self.ncBr5 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist3)
        self.ncBr5.weight[0] =0.0025 #changed from 0.002
        self.ncBr5.delay = 1
        self.nc2Br5 = h.NetCon(self._vecstim2, self.IN_1.syn_Brdist3)
        self.nc2Br5.weight[0] =0.0025 #changed from 0.002
        self.nc2Br5.delay = 1
        
        self.ncBr6 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist4)
        self.ncBr6.weight[0] =0.0025 #changed from 0.002
        self.ncBr6.delay = 1
        self.nc2Br6 = h.NetCon(self._vecstim2, self.IN_1.syn_Brdist4)
        self.nc2Br6.weight[0] =0.0025 #changed from 0.002
        self.nc2Br6.delay = 1
        
        # self.ncBr7 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist5)
        # self.ncBr7.weight[0] =0.0025 #changed from 0.002
        # self.ncBr7.delay = 1
        self.nc2Br7 = h.NetCon(self._vecstim2, self.IN_1.syn_Brdist5)
        self.nc2Br7.weight[0] =0.0025 #changed from 0.002
        self.nc2Br7.delay = 1

    def _create_cells(self, r):
        self.TCcells = []
        self.IN_1 = IN(self.gidStart+0,0,0,0,0,self.gidStart+1,self.gidStart+2,self.gidStart+3,self.gidStart+4)
        for i in range(self.gidStart+5,self._N+self.gidStart+5):    # only create the cells that exist on this host
            theta = i * 2 * h.PI / self._N
            self.TCcells.append(TC(i, h.cos(theta) * r, h.sin(theta) * r, 0, theta))
                      
    def _add_IN_syn(self):
        for target in self.TCcells:
            self.nc = h.NetCon(self.IN_1.model.soma[0](0.5)._ref_v, target.syn_INaxonal, sec=self.IN_1.model.soma[0])
            self.nc.weight[0] = 0.002 
            self.nc.delay = 0.97
            self.IN_1.ncs.append(self.nc)
            
        # for target in self.TCcells:
        #     if target._gid == 5 or target._gid == self.gidStart+5:
        #         self.nc2 = h.NetCon(self.IN_1.model.dend[20](1)._ref_v, target.syn_INdendritic, sec=self.IN_1.model.dend[20])
        #         self.nc2.weight[0] = 0.002 
        #         self.nc2.delay = 0.1
        #         self.nc2.threshold = -34 
        #         self.IN_1.ncs2.append(self.nc2)
        #     if target._gid == 6 or target._gid == self.gidStart+6:
        #         self.nc4 = h.NetCon(self.IN_1.model.dend[79](1)._ref_v, target.syn_INdendritic, sec=self.IN_1.model.dend[79])
        #         self.nc4.weight[0] = 0.002
        #         self.nc4.delay = 0.1
        #         self.nc4.threshold = -34 
        #         self.IN_1.ncs4.append(self.nc4)
        #     if target._gid == 7 or target._gid == self.gidStart+7:
        #         self.nc5 = h.NetCon(self.IN_1.model.dend[92](1)._ref_v, target.syn_INdendritic, sec=self.IN_1.model.dend[92])
        #         self.nc5.weight[0] = 0.002
        #         self.nc5.delay = 0.1
        #         self.nc5.threshold = -34 
        #         self.IN_1.ncs5.append(self.nc5)
        #     if target._gid == 8 or target._gid == self.gidStart+8:
        #         self.nc6 = h.NetCon(self.IN_1.model.dend[104](1)._ref_v, target.syn_INdendritic, sec=self.IN_1.model.dend[104])
        #         self.nc6.weight[0] = 0.002
        #         self.nc6.delay = 0.1
        #         self.nc6.threshold = -34 
        #         self.IN_1.ncs6.append(self.nc6)
            
    def _add_Br_Ps_syn(self):
        for target in self.TCcells:
            if target._gid == 5 or target._gid == self.gidStart+5:
                self.ncBrPs = h.NetCon(self._vecstim, target.syn_Br)
                self.ncBrPs.weight[0] = 0.024 #changed from 0.003
                self.ncBrPs.delay = 1
                self.ncBrsPs.append(self.ncBrPs)
            if target._gid == 6 or target._gid == self.gidStart+6:
                self.ncBrPs = h.NetCon(self._vecstim, target.syn_Br)
                self.ncBrPs.weight[0] = 0.024 #changed from 0.003
                self.ncBrPs.delay = 1
                self.ncBrsPs.append(self.ncBrPs)
                self.nc2BrPs = h.NetCon(self._vecstim2, target.syn_Br)
                self.nc2BrPs.weight[0] = 0.024 #changed from 0.003
                self.nc2BrPs.delay = 1
                self.nc2BrsPs.append(self.ncBrPs)
            if target._gid == 7 or target._gid == self.gidStart+7:
                self.ncBrPs = h.NetCon(self._vecstim, target.syn_Br)
                self.ncBrPs.weight[0] = 0.024 #changed from 0.003
                self.ncBrPs.delay = 1
                self.ncBrsPs.append(self.ncBrPs)
                self.nc2BrPs = h.NetCon(self._vecstim2, target.syn_Br)
                self.nc2BrPs.weight[0] = 0.024 #changed from 0.003
                self.nc2BrPs.delay = 1
                self.nc2BrsPs.append(self.ncBrPs)
            if target._gid == 8 or target._gid == self.gidStart+8:
                self.nc2BrPs = h.NetCon(self._vecstim2, target.syn_Br)
                self.nc2BrPs.weight[0] = 0.024 #changed from 0.003
                self.nc2BrPs.delay = 1
                self.nc2BrsPs.append(self.ncBrPs)