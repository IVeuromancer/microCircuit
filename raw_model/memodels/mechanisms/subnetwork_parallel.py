from neuron import h
from cell import IN
from cell import TC

h.nrnmpi_init()
pc = h.ParallelContext()

class Subnetwork:
    """
    A network of TC and IN cells where IN makes an
    axonal and dendritic synapse onto all other TCs and all cells 
    receive Br input
    
    """
    def __init__(self, gidStart=0, N=2, r=100, PS = None):
        """
        :param N: Number of cells. 2 by default
        :param r: radius of the network. 100um by default

        """
        self._N = N
        self.PS = PS #PS stands for poisson input
        self.gidStart = gidStart #this is to offset gids every time instantiate a subnetwork (e.g. in subnetwork1, IN has gids 0,1,2 and TCs have gids 3,4. Then subnetwork2 must start with gid 5)
        # Add simple Br spikes, these simple Br spikes are optional, I use vecstim for real simulation runs and turn off netstim before the run by making ._netstim.number = 0. 
        self._netstim = h.NetStim()
        # Br simple input parameters
        self._netstim.number = 10
        self._netstim.start = 1
        self._netstim.interval = 3
        self.ncBrs = []
        
        # Add Br spikes of nonstationary poisson process. from poisson.py
        self._vecstim = h.VecStim()
        self.train_vec = h.Vector(PS)
        self._vecstim.play(self.train_vec)
        self.ncBrsPs = []
        
        # Round-robin counting.
        # Each host has an id from 0 to pc.nhost() - 1.
        self.gidlist = list(range(pc.id()+gidStart, self._N+3, pc.nhost()))
        for gid in self.gidlist:
            pc.set_gid2node(gid+gidStart, pc.id())
        
        # create TCs and IN class objects and add synapses from IN to TC
        self._create_cells(r)
        
        # pc.multisplit(self.IN_1.model.dend[6](0.5), self.gidStart+10+1, 1)
        # pc.multisplit(self.IN_1.model.dend[16](0.5), self.gidStart+10+2, 1)

        # associate the cell with this host and gid
        pc.cell(self.IN_1._gid, self.IN_1._spike_detector)
        pc.cell(self.IN_1._gid2, self.IN_1._release_detector)
        pc.cell(self.IN_1._gid3, self.IN_1._release_detector2)
        for cell in self.TCcells:
            pc.cell(cell._gid, cell._spike_detector)
        
        # create NetCon from simple spikes and poisson process (external inputs) to TC
        self._add_IN_syn()
        self._add_Br_syn()
        self._add_Br_Ps_syn()
        
        # create NetCon from Br for IN targets (prox and dist denrites)
        self.ncBr2 = h.NetCon(self._netstim, self.IN_1.syn_Brprox)
        self.ncBr2.weight[0] =0.0008
        self.ncBr2.delay = 1

        self.ncBr3 = h.NetCon(self._netstim, self.IN_1.syn_Brdist)
        self.ncBr3.weight[0] =0.004
        self.ncBr3.delay = 1
        
        self.ncBrPs2 = h.NetCon(self._vecstim, self.IN_1.syn_Brprox)
        self.ncBrPs2.weight[0] =0.0008
        self.ncBrPs2.delay = 1

        self.ncBrPs3 = h.NetCon(self._vecstim, self.IN_1.syn_Brdist)
        self.ncBrPs3.weight[0] =0.004
        self.ncBrPs3.delay = 1

    def _create_cells(self, r):
        self.TCcells = []
        self.IN_1 = IN(self.gidStart+0,0,0,0,0,self.gidStart+1,self.gidStart+2)
        for i in range(self.gidStart+3,self._N+self.gidStart+3):    # only create the cells that exist on this host
            theta = i * 2 * h.PI / self._N
            self.TCcells.append(TC(i, h.cos(theta) * r, h.sin(theta) * r, 0, theta))
            
    def _add_IN_syn(self):
        for target in self.TCcells:
            self.nc = pc.gid_connect(self.IN_1._gid, target.syn_INaxonal)
            self.nc.weight[0] = 0.008
            self.nc.delay = 1
            self.nc.threshold = -15
            self.IN_1.ncs.append(self.nc)

        for target in self.TCcells:
            self.nc2 = pc.gid_connect(self.IN_1._gid2, target.syn_INdendritic)
            self.nc2.weight[0] = 0.008
            self.nc2.delay = 1
            self.nc2.threshold = -15 # threshold needs to be lower for distal dendritic release
            self.IN_1.ncs2.append(self.nc2)
            
    def _add_Br_syn(self): #optional simple Br spike input in case I need it
        for target in self.TCcells:
            self.ncBr = h.NetCon(self._netstim, target.syn_Br)
            self.ncBr.weight[0] = 0.04
            self.ncBr.delay = 1
            self.ncBrs.append(self.ncBr)
            
    def _add_Br_Ps_syn(self): #vecstim Br poisson spike input
        for target in self.TCcells:
            self.ncBrPs = h.NetCon(self._vecstim, target.syn_Br)
            self.ncBrPs.weight[0] = 0.04
            self.ncBrPs.delay = 1
            self.ncBrsPs.append(self.ncBrPs)