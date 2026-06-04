from neuron import h

class Cell:
    """
    Class objects for individual cells (either of type IN or TC)
    
    """
    def __init__(self, gid, x, y, z, theta, gid2 = None, gid3 = None):
        """
        :param gid: assignment of primary gid, every cell in the network has a unique gid
        :param x,y,z,theta: these position numbers are defined within the Subnetwork class, def _create_cells(). Position of the cells is not very important 
        :param gid2: assignment of second gid to IN type neurons only. dendrodendritic synapse to TC
        :param gid3: assignment of third gid to IN type neruons only. dendrodendritic synapse to other IN

        """
        self._gid = gid
        self._gid2 = gid2
        self._gid3 = gid3
        self._setup_cell()
        self.x = self.y = self.z = 0                     
        h.define_shape()
        self._rotate_z(theta)                                   
        self._set_position(x, y, z)  
        
        # record spikes, soma/prox-dend from all cells
        # ._spike_detector is identified by ._gid
        self._spike_detector = h.NetCon(self.model.soma[0](0.5)._ref_v, None, sec=self.model.soma[0])
        self._spike_detector.threshold = -15
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        self.soma_v = h.Vector().record(self.model.soma[0](0.5)._ref_v)
        self.dend_prox_v = h.Vector().record(self.model.dend[0](0.1)._ref_v)
    
        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
    def _set_position(self, x, y, z):
        for sec in self.model.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
    def _rotate_z(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.model.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))
    
class IN(Cell):
    name = 'Interneuron'
    def _setup_cell(self):
        h.load_file("../bAC_IN_legacy.hoc") 
        morphology_dir = "../morphologies_IN_bAC"
        morphology_name = "dend-jy180406_B_idC_axon-jy171019_B_10x_resta_idB_-_Clone_4.asc"
        self.model = h.bAC_IN_legacy(morphology_dir, morphology_name)
        
        # contains all spike detectors
        self.ncs = []

        # recording of release events in distal dendrites (there are dendro-dendritic synapses)
        # ._release_detector and ._release_detector2 are identified by ._gid2 and ._gid3
        self._release_detector = h.NetCon(self.model.dend[7](0.9)._ref_v, None, sec=self.model.dend[7])
        self._release_detector2 = h.NetCon(self.model.dend[17](0.9)._ref_v, None, sec = self.model.dend[17])
        self._release_detector.threshold = -15
        self._release_detector2.threshold = -15
        self.release_times = h.Vector()
        self.release_times2 = h.Vector()
        self._release_detector.record(self.release_times)
        self._release_detector2.record(self.release_times2)
        # contains all release detectors
        self.ncs2 = []
        self.dend_dist_v = h.Vector().record(self.model.dend[7](0.9)._ref_v)
        self.dend_dist2_v = h.Vector().record(self.model.dend[17](0.9)._ref_v)
        
        # add synapses received from brainstem Br
        self.syn_Brprox = h.Exp2Syn(self.model.dend[0](0.1))
        self.syn_Brprox.tau1 = 1.6
        self.syn_Brprox.tau2 = 3.6
        self.syn_Brprox.e = 10
        
        self.syn_Brdist = h.Exp2Syn(self.model.dend[7](0.9))
        self.syn_Brdist.tau1 = 0.3
        self.syn_Brdist.tau2 = 2
        self.syn_Brprox.e = 10
        
        # add current clamp to be able to inject current from all my cells, this is to depolarize the cells a bit
        self.stim = h.IClamp(self.model.soma[0](0.5))
        self.stim.dur = 0
        self.stim.delay = 0
        self.stim.amp = 0
        self.stim_current = h.Vector()
        self.stim_current.record(self.stim._ref_i)

class TC(Cell):
    name = 'TC'
    def _setup_cell(self):
        h.load_file("../dAD_ltb_legacy.hoc")
        morphology_dir = "../morphologies_TC_dAD_ltb"
        morphology_name = "dend-jy170925_B_idB_axon-AA0348_-_Scale_x1.000_y1.050_z1.000_-_Clone_2.asc"
        self.model = h.dAD_ltb_legacy(morphology_dir, morphology_name)
        
        # add synapses received from brainstem Br
        self.syn_Br = h.Exp2Syn(self.model.dend[0](0.1))
        self.syn_Br.tau1 = 0.2
        self.syn_Br.tau2 = 1.2
        self.syn_Br.e = 1
        
        # add synapses received from IN
        # INaxon to TC is .syn_INaxonal
        # INdendrite to TC is .syn_INdendritic (dendro-dendritic synapse)
        self.syn_INaxonal = h.Exp2Syn(self.model.dend[0](0.1))
        self.syn_INaxonal.tau1 = 0.7
        self.syn_INaxonal.tau2 = 4.2
        self.syn_INaxonal.e = -80
        
        self.syn_INdendritic = h.Exp2Syn(self.model.dend[0](0.1))
        self.syn_INdendritic.tau1 = 0.7
        self.syn_INdendritic.tau2 = 4.2
        self.syn_INdendritic.e = -80
        
        # add current clamp to be able to inject current from all my cells, this is to depolarize the cells a bit
        self.stim = h.IClamp(self.model.soma[0](0.5))
        self.stim.dur = 100
        self.stim.delay = 0
        self.stim.amp = 0.195
        self.stim_current = h.Vector()
        self.stim_current.record(self.stim._ref_i)
