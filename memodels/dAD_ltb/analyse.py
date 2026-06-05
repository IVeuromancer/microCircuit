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

import sys
import traceback
import bluepyopt
import os
import collections

import matplotlib
matplotlib.use('Agg', warn=True)
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

def str2bool(v):
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "0"):
        return False
    elif v.lower() in ("none"):
        return None
    else:
        return v

import tools
import json

#main_path = os.getcwd()
#main_path = os.path.dirname(__file__)
main_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(main_path)

import argparse
parser = argparse.ArgumentParser(description='post')

parser.add_argument('--rundir', required=False, default=".")
parser.add_argument('--githash', type=str, required=False, default="None")
parser.add_argument('--seed', type=str, required=False, default='1')
parser.add_argument('--rank', type=int, required=False, default=0)
parser.add_argument('--branch', type=str, required=False, default="None")

parser.add_argument('--nosimulate', action="store_false")
parser.add_argument('--test', required=False, default="None")

parser.add_argument('--stochdet', type=str, required=False, default="None")
parser.add_argument('--usethreshold', action="store_true")
parser.add_argument('--mmtest', action="store_true",
                help='Run different recipe mm_test_recipe as give in recipes.json')

parser.add_argument('--hoc', required=False, default=False,
                help='Supply a hoc file to be executed')
parser.add_argument('--oldhoc', required=False, default=False,
                help='Supply an old style hoc file to be executed')
parser.add_argument('--nrnivmodl', action="store_true",
                help='Recompile mod files')
parser.add_argument('--altmorph', type=str, required=False, default="None")


args, unknown = parser.parse_known_args()

exec('seeds=['+args.seed +']')
githash = str2bool(args.githash)
rundir = args.rundir
rank = args.rank
branch = args.branch
simulate = args.nosimulate
etypetest = str2bool(args.test)
hoc = args.hoc
oldhoc = args.oldhoc
nrnivmodl = args.nrnivmodl
stochdet = str2bool(args.stochdet)
usethreshold = args.usethreshold
altmorph = str2bool(args.altmorph)
mmtest = args.mmtest


def main():
    """Main"""

    global usethreshold

    figs=collections.OrderedDict()
    grouping = ['etype', 'githash', 'seed', 'rank', 'stochdet', 'altmorph']

    
    for seed in seeds:

        analyse = tools.Analyse(rundir=rundir, githash=githash, seed=seed, rank=rank,
                            etypetest=etypetest, hoc=hoc, oldhoc=oldhoc,
                            nrnivmodl=nrnivmodl, grouping=grouping, main_path=main_path,
                            altmorph=altmorph)

        reportname, label = analyse.get_name()
        reportname = reportname.split('seed')[0]

        if mmtest or (branch is not "None"):
            recipe = json.load(open('config/recipes/recipes.json'))
            if 'mm_test_recipe' in recipe[analyse.etype]:
                analyse.etype = recipe[analyse.etype]['mm_test_recipe']
                usethreshold = False

        threshold_fit = ("_legacy" not in reportname) and ("_combined" not in reportname)

        do_usethreshold = usethreshold and (not threshold_fit)
        emodel = analyse.etype[:]

        if do_usethreshold:
            split_sim = 2
        else:
            split_sim = 1

        color = 'b'
        analyse.plot_evolution(figs=figs, color=color)
        analyse.do_model_export()
        analyse.plot_diversity(figs=figs, color=color)
        analyse.save_pdf(reportname, figs)

        analyse.sim_plot_obj(stochdet=None, usethreshold=False,
                    figs=figs, color=color, split_sim=split_sim)
        analyse.save_pdf(reportname, figs)

        if branch is not "None":
            analyse.add_to_final(branch, emodel)

        if do_usethreshold:
            color='r'
            analyse.sim_plot_obj(stochdet=None, usethreshold=True,
                        figs=figs, color=color, split_sim=split_sim)
            analyse.save_pdf(reportname, figs)


        if ((stochdet is False)
            and (('IR' in reportname) or ('STUT' in reportname))
            and ('_det' not in reportname)):

            color='b'
            analyse.sim_plot_obj(stochdet=False, usethreshold=False,
                        figs=figs, color=color, split_sim=split_sim)
            analyse.save_pdf(reportname, figs)

            if do_usethreshold:
                color='r'
                analyse.sim_plot_obj(stochdet=False, usethreshold=True,
                            figs=figs, color=color, split_sim=split_sim)
            analyse.save_pdf(reportname, figs)


    plt.show()


if __name__ == '__main__':
    try:
        main()
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
