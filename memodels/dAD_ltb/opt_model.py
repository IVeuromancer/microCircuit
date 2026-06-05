#!/usr/bin/env python

"""Run simple cell optimisation"""

"""
Copyright (c) 2016, EPFL/Blue Brain Project

 This file is part of BluePyOpt <https://github.com/BlueBrain/BluePyOpt>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

# pylint: disable=R0914

import sys
import traceback
import bluepyopt
import efel
import os

# TODO rename 'score' into 'objective'
# TODO add functionality to read settings of every object from config format

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

import setup  # pylint: disable=W0403
import json
#import yaml
import collections
import subprocess

if os.environ.get("USE_NEXUS_FOR_THALAMUS") is not None:
    from input_from_nexus import recipes_from_nexus, register_current_run
    recipes = json.load(open(recipes_from_nexus()))
    etypes = recipes.keys()
else:
    recipes = json.load(open('config/recipes/recipes.json'))
    etypes = recipes.keys()

def makedirs(filename): # also accepts filename
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

def str2bool(v):
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "0"):
        return False
    elif v.lower() in ("none"):
        return None
    else:
        return v






import argparse
parser = argparse.ArgumentParser(description='cell')
parser.add_argument('--start', action="store_true")
parser.add_argument('--resume', action="store_true", default=False)

parser.add_argument('--etype', required=False, choices=etypes, default="cADpyr_L5TPC")
parser.add_argument('--offspring_size', type=int, required=False, default=2,
                    help='number of individuals in offspring')
parser.add_argument('--max_ngen', type=int, required=False, default=100,
                    help='maximum number of generations')
parser.add_argument('--seed', type=int, required=False, default=1,
                    help='seed for optimization')

parser.add_argument('--githash', type=str, required=False, default="None")
parser.add_argument('--runlog', required=False, default="./runlog.json")
parser.add_argument('--slurmid', type=int, required=False, default=0)

parser.add_argument('--stochdet', type=str, required=False, default="None")
parser.add_argument('--usethreshold', action="store_true")

args, unknown = parser.parse_known_args()

seed = args.seed
etype = args.etype
offspring_size = args.offspring_size
max_ngen = args.max_ngen
githash = str2bool(args.githash)
runlogpath = args.runlog
slurmid = args.slurmid
stochdet = str2bool(args.stochdet)
usethreshold = args.usethreshold

if os.environ.get("USE_NEXUS_FOR_THALAMUS") is not None:
    register_current_run(githash, seed)

if args.resume:

    metafile = './meta/%s_%d.yaml' % (githash, seed)
    if os.path.isfile(metafile):
        meta = yaml.load(open(metafile, 'r'))
        etype = meta['etype']
        offspring_size = meta['offspring_size']
        notes = ""
    else:
        runlog = json.load(open(runlogpath, 'r'))
        meta = runlog[str(githash)][str(seed)]
        etype = meta['etype']
        offspring_size = meta['offspring_size']
        notes = meta['notes']

    #if 'githash' in meta:
    #    subprocess.check_output(["git", "checkout", meta['githash']])

elif args.start:

    import setup.tools as tools

    # Raises an IOError in 5 seconds if unable to acquire the lock.
    with tools.SimpleFlock(runlogpath, timeout=5):
        try:
            with open(runlogpath) as f:
                runlog = json.load(f,
                    object_pairs_hook=collections.OrderedDict)
        except:
            runlog = collections.OrderedDict()

        if str(githash) not in runlog:
            runlog[str(githash)] = collections.OrderedDict()
        # else:
        #     runlog[str(githash)] = collections.OrderedDict(sorted(runlog[str(githash)].items()))

        runlog[str(githash)][str(seed)] = collections.OrderedDict([
                                        ('githash',githash),
                                        ('slurmid',slurmid),
                                        ('etype',etype),
                                        ('offspring_size',offspring_size),
                                        ('notes',''),
                                        ('bluepyopt_version',bluepyopt._version.get_versions()),
                                        ('efel_version',efel._version.get_versions()),
                                    ])

        s = json.dumps(runlog, indent=4)
        s = tools.collapse_json(s, indent=8)

        with open(runlogpath, "w") as f:
            f.write(s)

        with open(os.path.basename(runlogpath), "w") as f:
            f.write(s)

evaluator = setup.evaluator.create(etype=etype,
                                   stochkv_det=stochdet,
                                   usethreshold=usethreshold,
                                   runopt=True)

checkpoints_dir = './checkpoints/run.%s' % githash
cp_filename = os.path.join(
    checkpoints_dir, 'checkpoint_%d.pkl' % seed)

if os.getenv('USEIPYP') == '1':
    from ipyparallel import Client
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
    lview = rc.load_balanced_view()

    map_function = lview.map_sync
    live_plot = False
else:
    map_function = None
    live_plot = True

import setup.protocols


opt = bluepyopt.optimisations.DEAPOptimisation(
    evaluator=evaluator,
    map_function=map_function,
    seed=seed,
    eta=10., mutpb=1.0, cxpb=1.0)


def main():
    """Main"""
    print args
    if args.start or args.resume:
        logger.debug('Doing start or continue')

        #makedirs(checkpoints_dir)
        opt.run(max_ngen=max_ngen,
                offspring_size=offspring_size,
                continue_cp=args.resume,
                cp_filename=cp_filename)


if __name__ == '__main__':
    try:
        main()
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
