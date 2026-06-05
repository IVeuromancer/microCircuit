
from bluepyopt.ephys.models import CellModel

import bglibpy
from bglibpy.importer import neuron

import collections
import logging

logger = logging.getLogger(__name__)
import sys

class HocOldModel(CellModel):

    '''Wrapper class for an old hoc template so it can be used by BluePyOpt'''

    def __init__(self, name, morphology_path, hoc_path):
        """Constructor

        Args:
            name(str): name of this object
            sim(NrnSimulator): simulator in which to instatiate hoc_path
            morphology_path(str path): path to morphology that can be loaded by
                                       Neuron
            hoc_path(str path): path to .hoc file that will be used
        """
        super(HocOldModel, self).__init__(name,
                                           morph=None,
                                           mechs=[],
                                           params=[])
        self.hoc_path = hoc_path
        self.morphology_path = morphology_path
        self.cell = None
        self.icell = None

    def params_by_names(self, param_names):
        pass

    def freeze(self, param_dict):
        pass

    def unfreeze(self, param_names):
        pass

    def instantiate(self, sim=None):
        """Instantiate model in bglibpy"""

        self.cell = bglibpy.Cell(self.hoc_path, self.morphology_path)

        try:
            self.cell.re_init_rng(use_random123_stochkv=True)
        except AttributeError:
            sys.exc_clear()

        self.name = self.cell.cellname
        self.icell = self.cell.cell.getCell()

    def destroy(self, sim=None):
        self.cell = None
        self.icell = None

    def check_nonfrozen_params(self, param_names):
        pass

    def __str__(self):
        """Return string representation"""
        return ('%s: %s of %s(%s)' %
                (self.__class__, self.name, self.hoc_path,
                 self.morphology_path, ))
