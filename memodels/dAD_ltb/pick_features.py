#from neo import io
import numpy
import sys
import os
import fnmatch

try:
    import cPickle as pickle
except:
    import pickle

import bluepyefe as bpefe

import gzip
import json
from collections import OrderedDict
from copy import deepcopy

# Don't overwrite atm
"""
import shutil
try:
    shutil.rmtree('./config/features')
    shutil.rmtree('./config/protocols')
except:
    pass
"""
delay = 0

mainpath = "/gpfs/bbp.cscs.ch/project/proj55/singlecell/singlecell-features/mouse"

RMP = OrderedDict()
RMP["type"] = "StepProtocol"
RMP["stimuli"] = {}
RMP["stimuli"]["step"] = OrderedDict([ ("delay", 300), ("amp", 0), ("duration", 900.0), ("totduration", 1200.0) ])

ThProt = OrderedDict()
ThProt["Main"] = OrderedDict()
ThProt["Main"]["type"] = "RatSSCxMainProtocol"
ThProt["Main"]["rmp_score_threshold"] = 3
ThProt["Main"]["rin_score_threshold"] = 3
ThProt["Main"]["other_protocols"] = []
ThProt["Main"]["pre_protocols"] = []
ThProt["Main"]["preprot_score_threshold"] = 3

for suffix in {"hyp", "dep"}:
    ThProt["RinHoldcurrent_{}".format(suffix)] = OrderedDict()
    ThProt["RinHoldcurrent_{}".format(suffix)]["type"] = "RatSSCxRinHoldcurrentProtocol"
    ThProt["RinHoldcurrent_{}".format(suffix)]["holdi_precision"] = 0.1
    ThProt["RinHoldcurrent_{}".format(suffix)]["holdi_max_depth"] = 14
    ThProt["RinHoldcurrent_{}".format(suffix)]["holdi_estimate_multiplier"] = 1
    ThProt["ThresholdDetection_{}".format(suffix)] = OrderedDict()
    ThProt["ThresholdDetection_{}".format(suffix)]["type"] = "RatSSCxThresholdDetectionProtocol"
    ThProt["ThresholdDetection_{}".format(suffix)]["step_template"] = OrderedDict()
    ThProt["ThresholdDetection_{}".format(suffix)]["step_template"]["type"] = "StepProtocol"
    ThProt["ThresholdDetection_{}".format(suffix)]["step_template"]["stimuli"] = OrderedDict()
    ThProt["ThresholdDetection_{}".format(suffix)]["step_template"]["stimuli"]["step"] = OrderedDict([ ("delay", 800.0), ("amp", None), ("duration", 1800.0), ("totduration", 1800.0) ])
    ThProt["ThresholdDetection_{}".format(suffix)]["step_template"]["stimuli"]["holding"] = OrderedDict([ ("delay", 0.0), ("amp", None), ("duration", 1800.0), ("totduration", 1800.0) ])

ThProt_val = deepcopy(ThProt)
ramp_delay = 800
for suffix in {"_hyp", ""}:
    ThProt_val["Ramp"+suffix] = OrderedDict()
    ThProt_val["Ramp"+suffix]["type"] = "RampThresholdProtocol"
    ThProt_val["Ramp"+suffix]["stimuli"] = {}
    dur = 4000 if (suffix == "") else 1250
    ThProt_val["Ramp"+suffix]["stimuli"]["ramp"] = OrderedDict({"delay": 800.0, "amp": None, "thresh_perc_start": 0.0, "thresh_perc_end": 300.0, "duration": dur, "totduration": dur+ramp_delay})

getfeat = {
    'bAC_IN':[
        {   'source': 'IN/bAC', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_200','IDRest_250'],#'IDRest_300'], # Only 2 cells at 300
            'rename': ['Step_150_hyp','Step_200_hyp','Step_250_hyp'],#'Step_300_hyp'],
            'features':[
                'voltage_after_stim', # v
                'AP_amplitude', 'APlast_amp', 'AHP_depth', # AP
                'mean_frequency', 'inv_time_to_first_spike', 'time_to_last_spike', 'inv_first_ISI', 'inv_second_ISI', # freq
                'inv_third_ISI', 'inv_fourth_ISI', 'inv_fifth_ISI', 'inv_last_ISI', # freq
            ],
        },
        {   'source': 'IN/bAC', 'location':['soma.v'],
            'names':['IV_-40'],
            'rename':['Rin_hyp'],
            'features':['ohmic_input_resistance_vb_ssse', 'voltage_base'],
        },
        {   'source': 'IN/bAC', 'location':['soma.v'],
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'IN/bAC', 'location':['soma.v'],
            'names':['IV_-120'],
            'rename':['IV_-120_hyp'],
            'features':['voltage_deflection', 'voltage_deflection_begin'],
        }
    ],
    'cNAD_noscltb':[
        {   'source': 'Rt/cNAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'],
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                #'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                "initburst_sahp",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':['IV_-120'], # -20 for Rin
            'features':['voltage_deflection'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base','ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ],
    'cAD_noscltb':[
        {   'source': 'Rt/cAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'],
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                #'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                "initburst_sahp",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':['IV_-120'], # -20 for Rin
            'features':['voltage_deflection'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base','ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'Rt/noscltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ],
    'cAD_oscltb':[
        {   'source': 'Rt/cAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'],
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                #'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'Rt/oscltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'burst_number',
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                "initburst_sahp",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':['IV_-120'], # -20 for Rin
            'features':['voltage_deflection'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'Rt/oscltb', 'location':['soma.v'],
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/oscltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/cAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base','ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'Rt/oscltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ],
    'cNAD_nltb':[
        {   'source': 'Rt/cNAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'],
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                #'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'Rt/nltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                "initburst_sahp",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':['IV_-120'], # -20 for Rin
            'features':['voltage_deflection'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'Rt/nltb', 'location':['soma.v'],
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/nltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'Rt/cNAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'Rt/nltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ],
    'dNAD_ltb':[
        {   'source': 'TC/dNAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'TC/dNAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'], # No IDRest 250
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'TC/dNAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                #"initburst_sahp_vb",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'TC/dNAD', 'location':['soma.v'],
            'names':['IV_-140'], # -20 for Rin
            'features':['sag_amplitude'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'TC/ltb', 'location':['soma.v'], # There are more recordings in bursting mode, TODO: RMP should be average from all recordings
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'TC/dNAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'TC/dNAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ],
    'dAD_ltb':[
        {   'source': 'TC/dAD', 'location':['soma.v'], # for hypamp and threshold, dummy
            'names':[]
        },
        {   'source': 'TC/dAD', 'location':['soma.v'],
            'names': ['IDRest_150','IDRest_250'], # No IDRest 250
            'rename': ['Step_150', 'Step_250'],
            'features':[
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before
                'Spikecount', 
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI',  'inv_last_ISI',# freq, removed time_to_first_spike, removed mean_freq
                'adaptation_index2' # freq
            ],
        },
        {   'source': 'TC/dAD', 'location':['soma.v'],
            'names': ['IDRest_200'], # No IDRest 250
            'rename': ['Step_200'],
            'features':[
                'voltage_base',
                'AP_amplitude', 'AHP_depth', 'AP_duration_half_width', # 'AHP_depth_abs' before                  
                'Spikecount',
                'time_to_first_spike','inv_first_ISI', 
                'inv_second_ISI', # freq, removed time_to_first_spike, removed mean_freq
                'inv_last_ISI', # freq                        
                'adaptation_index2'
            ],
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names': ['IDRest_200'],
            'rename': ['Step_200_hyp'],
            'features':[
                'voltage_base',  # v
                #"AP2_AP1_peak_diff", # seems related to AHP after burst
                #'AP1_amp', 'AHP1_depth_from_peak', 'AHP2_depth_from_peak', # AP, Removed: 'amp_drop_first_second', 
                #'AP1_width', 'AP2_width',
                'AP1_amp', 'AP2_amp',
                #'AP_width', # width of LTS
                #"initburst_sahp_vb",
                'time_to_first_spike', 'Spikecount' ,  # freq
                'inv_first_ISI', 'inv_second_ISI', 'inv_last_ISI', # freq
                'voltage_after_stim' # v
            ],

        },
        {   'source': 'TC/dAD', 'location':['soma.v'],
            'names':['IV_-140'], # -20 for Rin
            'features':['sag_amplitude'
                         ], # Before: time_constant, 'voltage_deflection', 'voltage_deflection_begin', 'ohmic_input_resistance_vb_ssse'
        },
        {   'source': 'TC/ltb', 'location':['soma.v'], # There are more recordings in bursting mode, TODO: RMP should be average from all recordings
            'names':["SponNoHold30_all"],
            'rename':["RMP"],
            'features':['voltage_base', 'Spikecount'],
            'featrename':['steady_state_voltage_stimend', 'Spikecount_stimint'],
            'protocol': RMP
        },
        {   'source': 'TC/dAD', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_dep"],
            'features':['Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names':["SponHold30_all"], # TODO: check which is the best protocol for this
            'rename':["hold_hyp"],
            'features':[ 'Spikecount'], # 'voltage_base',
            'featrename':['Spikecount_stimint'], # 'steady_state_voltage_stimend',
        },
        {   'source': 'TC/dAD', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_dep'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        },
        {   'source': 'TC/ltb', 'location':['soma.v'],
            'names':['IV_-40'], # Before -40 for Rin, -140 for sag
            'rename':['Rin_hyp'],
            'features':['voltage_base', 'ohmic_input_resistance_vb_ssse' ],
        }

    ]

}

for etypename, etype in getfeat.iteritems():

    etype_features = OrderedDict()
    etype_protocols = OrderedDict()

    ThProt_ = deepcopy(ThProt)
    if '_small' in etypename:
        ThProt_["Main"]["rmp_score_threshold"] = 10
        ThProt_["Main"]["rin_score_threshold"] = 10
        ThProt_["RinHoldcurrent"]["holdi_precision"] = 0.1
        ThProt_["RinHoldcurrent"]["holdi_max_depth"] = 10

    etype_protocols_new = deepcopy(ThProt_)
    
    etype_protocols_new_val = deepcopy(ThProt_val)
    etype_protocols_new_val["Main"]["other_protocols"].append("Ramp")
    etype_protocols_new_val["Main"]["other_protocols"].append("Ramp_hyp")

    etype_protocols_combined = deepcopy(ThProt_)

    for group in etype:
        allfeatures = json.load(open(os.path.join(mainpath, group['source'], 'features.json')),
                        object_pairs_hook=OrderedDict)
        allprotocols = json.load(open(os.path.join(mainpath, group['source'], 'protocols.json')),
                        object_pairs_hook=OrderedDict)

        for iname, oldname in enumerate(group['names']):

            # protocols
            protocol = allprotocols[oldname]

            if 'rename' in group:
                name = group['rename'][iname]
            else:
                name = oldname

            if 'preprot_score_threshold' in group:
                etype_protocols_new["Main"]["preprot_score_threshold"] = group["preprot_score_threshold"]
                etype_protocols_combined["Main"]["preprot_score_threshold"] = group["preprot_score_threshold"]

            if 'location' in group:
                # protocols
                old_delay = protocol['stimuli']['step']['delay']
                protocol['stimuli']['step']['delay'] = old_delay

                #totduration = protocol['stimuli']['step']['totduration']
                totduration = old_delay + 250.0 + protocol['stimuli']['step']['duration']
                protocol['stimuli']['step']['totduration'] = totduration
                protocol['stimuli']['holding']['duration'] = totduration
                protocol['stimuli']['holding']['totduration'] = totduration
                #protocol['stimuli']['holding']['delay'] = totduration

                if 'protoptions' in group:
                    for opt, val in group['protoptions'].iteritems():
                        protocol['stimuli']['step'][opt] = val

                if 'stage' in group:
                    protocol['stage'] = group['stage'] # define fitting stage

                if 'protocol' in group:
                    etype_protocols[name] = group['protocol']
                    etype_protocols_new[name] = group['protocol']
                    etype_protocols_combined[name] = group['protocol']
                    if name in ["Rin_hyp", "Rin_dep", "RMP"]:
                        etype_protocols_new_val[name] = group['protocol']

                else:
                    etype_protocols[name] = deepcopy(protocol)

                    if name not in ['Rin_dep', 'Rin_hyp']:
                        if 'preprot' in group:
                            etype_protocols_new["Main"]["pre_protocols"].append(name)
                            etype_protocols_combined["Main"]["pre_protocols"].append(name)
                        else:
                            etype_protocols_new["Main"]["other_protocols"].append(name)
                            etype_protocols_combined["Main"]["other_protocols"].append(name)
                    else:
                        protocol['stimuli']['holding']['amp'] = None

                    etype_protocols_combined[name] = deepcopy(protocol)

                    # convert to threshold protocols
                    if (name is not 'bAP') and (name is not "Rin_dep" ) and (name is not "Rin_hyp" ):
                        # remove holding
                        del protocol['stimuli']['holding']
                        protocol['type'] = "StepThresholdProtocol"
                        protocol['stimuli']['step']['amp'] = None

                    if "threshold" in protocol['stimuli']['step']:
                        th = protocol['stimuli']['step']['threshold']
                        del protocol['stimuli']['step']['threshold']
                        protocol['stimuli']['step']['thresh_perc'] = th

                    etype_protocols_new[name] = deepcopy(protocol)
                    if name in ["Rin_hyp", "Rin_dep", "RMP"]:
                        etype_protocols_new_val[name] = deepcopy(protocol)
 
                # features
                etype_features[name] = OrderedDict()
                #for location, features in allfeatures[oldname].iteritems():
                for location in group['location']:
                    features = allfeatures[oldname][location]
                    order = []
                    subfeatures = []
                    for feature in features:

                        if feature['feature'] in group['features']:
                            if 'featrename' in group:
                                position = group['features'].index(feature['feature'])
                                feature['feature'] = group['featrename'][position]

                            if 'weights' in group:
                                if feature['feature'] in group['weights']:
                                    feature['weight'] = group['weights'][feature['feature']]

                            if 'stage' in group:
                                feature['stage'] = group['stage'] # define fitting stage

                            subfeatures.append(feature)

                    if 'featrename' in group:
                        order = group['featrename']
                    else:
                        order = group['features']
                    subfeatures = sorted(subfeatures, key=lambda x: order.index(x['feature']))

                    if 'extra' in group:
                        subfeatures += group['extra']

                    etype_features[name][location] = subfeatures

            else:
                if 'preprot' in group:
                    etype_protocols_new["Main"]["pre_protocols"].append(name)
                    etype_protocols_combined["Main"]["pre_protocols"].append(name)
                else:
                    etype_protocols_new["Main"]["other_protocols"].append(name)
                    etype_protocols_combined["Main"]["other_protocols"].append(name)

                etype_protocols[name] = protocol
                etype_protocols_new[name] = protocol
                etype_protocols_combined[name] = protocol
                etype_features[name] = allfeatures[oldname]

    # add hypamp and threshold features
    allhypth = json.load(open(os.path.join(mainpath, etype[0]['source'], 'hypamp_threshold.json')),
                    object_pairs_hook=OrderedDict)

    etype_features_new = deepcopy(etype_features)

    # if '_small' not in etypename:
    #     etype_features_new['RinHoldCurrent']= {'soma.v':[{'feature':'bpo_holding_current', 'val':allhypth["all"]["hypamp"] }]}
    #     etype_features_new['Threshold']= {'soma.v':[{'feature':'bpo_threshold_current', 'val':allhypth["all"]["threshold"] }]}

    if ('_abs' not in etypename) and ('_small' not in etypename):
        # These protocols are needed only in legacy
        # del etype_features_new["hold_dep"]
        # del etype_features_new["hold_hyp"]
        # del etype_protocols_new["hold_dep"]
        # del etype_protocols_new["hold_hyp"]
        # This protocol needed only if not legacy
        if "ltb" in etypename:
            del etype_features["Rin_dep"]
            del etype_protocols["Rin_dep"]

        # Overwrite spikecount in RMP protocol
        for feat_dict in [etype_features, etype_features_new]:
            for feature in feat_dict["RMP"]["soma.v"]:
                if "Spikecount_stimint" in feature.values() or "Spikecount" in feature.values() :
                    feature["val"] = [0.0, 0.001]
                    #feature["threshold"] = -55
                    #feature["max_score"] = 1000

        if "hold_dep" in etype_protocols.keys():
            for prot_dict in [etype_protocols, etype_protocols_new]:
                # Overwrite durations of these protocols
                prot_dict["hold_dep"]["stimuli"]["step"]["duration"] = 900.0
                prot_dict["hold_dep"]["stimuli"]["step"]["delay"] = 300.0
                prot_dict["hold_dep"]["stimuli"]["step"]["thresh_perc"] = 0.0
                prot_dict["hold_dep"]["stimuli"]["step"]["amp"] = 0.0
                prot_dict["hold_dep"]["stimuli"]["step"]["totduration"] = 1200.0
                try:
                    prot_dict["hold_dep"]["stimuli"]["holding"]["duration"] = 1200.0
                    prot_dict["hold_dep"]["stimuli"]["holding"]["totduration"] = 1200.0
                    prot_dict["hold_dep"]["stimuli"]["holding"]["amp"] = etype_protocols["Step_250"]["stimuli"]["holding"]["amp"]
                except KeyError:
                    pass
            # Overwrite spikecount in hold_dep protocol
            for feat_dict in [etype_features, etype_features_new]:
                for feature in feat_dict["hold_dep"]["soma.v"]:
                    if "Spikecount_stimint" in feature.values() or "Spikecount" in feature.values() :
                        feature["val"] = [0.0, 0.001]
                        feature["threshold"] = -55
                        feature["max_score"] = 1000
        if "hold_hyp" in etype_protocols.keys():
            for prot_dict in [etype_protocols, etype_protocols_new]:
                # Overwrite durations of these protocols
                prot_dict["hold_hyp"]["stimuli"]["step"]["duration"] = 900.0
                prot_dict["hold_hyp"]["stimuli"]["step"]["delay"] = 300.0
                prot_dict["hold_hyp"]["stimuli"]["step"]["thresh_perc"] = 0.0
                prot_dict["hold_hyp"]["stimuli"]["step"]["amp"] = 0.0
                prot_dict["hold_hyp"]["stimuli"]["step"]["totduration"] = 1200.0
                try:
                    prot_dict["hold_hyp"]["stimuli"]["holding"]["duration"] = 1200.0
                    prot_dict["hold_hyp"]["stimuli"]["holding"]["totduration"] = 1200.0
                    prot_dict["hold_hyp"]["stimuli"]["holding"]["amp"] = etype_protocols["Step_200_hyp"]["stimuli"]["holding"]["amp"]
                except KeyError:
                    pass
            # Overwrite spikecount in hold_dep protocol
            for feat_dict in [etype_features, etype_features_new]:
                for feature in feat_dict["hold_hyp"]["soma.v"]:
                    if "Spikecount_stimint" in feature.values() or "Spikecount" in feature.values() :
                        feature["val"] = [0.0, 0.001]
                        feature["threshold"] = -55
                        feature["max_score"] = 1000

        etype_features_path = os.path.join('./config/features/', etypename + '.json')
        bpefe.makedirs(etype_features_path)
        s = json.dumps(etype_features_new, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_features_path, "w") as f:
            f.write(s)

        etype_features_path = os.path.join('./config/features/', etypename + '_val.json')
        bpefe.makedirs(etype_features_path)
        etype_features_new_val = {prot:feat for (prot,feat) in etype_features_new.items()
                                  if "Rin" in prot or "RMP" in prot or "Threshold" in prot}
        # Overwrite feat for hyp hold search for validation with ramp
        for feature in etype_features_new_val["Rin_hyp"]["soma.v"]:
            if 'voltage_base' in feature.values():
                feature["val"][0], feature["n"] = -95.0, None

        etype_protocol_path = os.path.join('./config/protocols/', etypename + '.json')
        bpefe.makedirs(etype_protocol_path)
        s = json.dumps(etype_protocols_new, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_protocol_path, "w") as f:
            f.write(s)

        etype_protocol_path = os.path.join('./config/protocols/', etypename + '_val.json')
        bpefe.makedirs(etype_protocol_path)
        s = json.dumps(etype_protocols_new_val, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_protocol_path, "w") as f:
            f.write(s)

        etype_protocol_path = os.path.join('./config/protocols/', etypename + '_combined.json')
        bpefe.makedirs(etype_protocol_path)
        s = json.dumps(etype_protocols_combined, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_protocol_path, "w") as f:
            f.write(s)


        etype_features_path = os.path.join('./config/features/', etypename + '_legacy.json')
        bpefe.makedirs(etype_features_path)


        s = json.dumps(etype_features, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_features_path, "w") as f:
            f.write(s)

        etype_protocol_path = os.path.join('./config/protocols/', etypename + '_legacy.json')


        # Overwrite durations of these protocols
        try:
            etype_protocols["hold_dep"]["stimuli"]["step"]["duration"] = 900.0
            etype_protocols["hold_dep"]["stimuli"]["step"]["totduration"] = 1200.0
            etype_protocols["hold_dep"]["stimuli"]["holding"]["duration"] = 1200.0
            etype_protocols["hold_dep"]["stimuli"]["holding"]["totduration"] = 1200.0
            etype_protocols["hold_hyp"]["stimuli"]["step"]["duration"] = 900.0
            etype_protocols["hold_hyp"]["stimuli"]["step"]["totduration"] = 1200.0
            etype_protocols["hold_hyp"]["stimuli"]["holding"]["duration"] = 1200.0
            etype_protocols["hold_hyp"]["stimuli"]["holding"]["totduration"] = 1200.0
        except KeyError:
            pass

        bpefe.makedirs(etype_protocol_path)
        s = json.dumps(etype_protocols, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_protocol_path, "w") as f:
            f.write(s)

    else:
        etype_features_path = os.path.join('./config/features/', etypename + '.json')
        bpefe.makedirs(etype_features_path)
        s = json.dumps(etype_features, indent=2)
        s = bpefe.collapse_json(s, indent=6)
        with open(etype_features_path, "w") as f:
            f.write(s)

        etype_protocol_path = os.path.join('./config/protocols/', etypename + '.json')
        bpefe.makedirs(etype_protocol_path)
        s = json.dumps(etype_protocols, indent=2)
        s = bpefe.collapse_json(s, indent=6)

        with open(etype_protocol_path, "w") as f:
            f.write(s)
