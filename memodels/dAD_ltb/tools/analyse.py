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
import bluepyopt.ephys as ephys

import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import logging
logger = logging.getLogger(__name__)

import json
json.encoder.FLOAT_REPR = lambda x: format(x, '.17g')

#import yaml
from collections import OrderedDict
import subprocess
import pickle
import numpy
import time
import sh
import shutil
import plot
import plottools as pt


def makedirs(filename): # also accepts filename
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


class Analyse(object):

    # build cell
    def __init__(self, rundir, githash, seed=1, rank=0,
            etype=None, # use this to evaluate model
            etypetest=None, # also a test model can be given to evaluate against the model given in etype
            hoc=False, oldhoc=False, nrnivmodl=False,
            main_path=None,
            recipes_path='config/recipes/recipes.json',
            grouping=['etype', 'githash', 'seed', 'rank', 'altmorph'],
            figpath='figures',
            altmorph=None):

        self.githash = githash
        self.seed = seed
        self.rank = rank
        self.etypetest = etypetest
        self.etype = etype

        # options that can be used later:
        self.usethreshold = 'False'
        self.stochdet = None

        self.grouping = grouping
        self.notes = ""
        self.parameters = False

        self.hoc = hoc
        self.oldhoc = oldhoc

        self.altmorph = altmorph

        # switch to run dir and load modules
        path = os.path.join(main_path, rundir)
        sys.path.insert(0, path)
        os.chdir(path)
        import setup
        from setup.tools import collapse_json
        self.collapse_json = collapse_json
        self.setup = setup

        import inspect
        argspec = inspect.getargspec(self.setup.evaluator.create)
        self.has_etest = ('etypetest' in argspec.args)

        logger.info('Loading modules from %s', setup.__file__)
        #self.currentdir = currentdir

        if nrnivmodl:
            try:
               shutil.rmtree("x86_64")
            except:
               pass
            sh.nrnivmodl("mechanisms/")

        self.path_final = os.path.join(main_path, "final.json")
        self.main_path = main_path

        self.figpath = os.path.join(self.main_path, figpath)
        makedir(self.figpath)

        # check if etype is given to evaluate cell
        if self.etype is None:

            # try global fist
            runlogpath = os.path.join(main_path, "runlog.json")
            runlog = json.load(open(runlogpath, 'r'))

            # try local one if nothing in global
            if githash is not None:
                if str(githash) not in runlog:
                    runlogpath = "runlog.json"
                    runlog = json.load(open(runlogpath, 'r'))

        self.recipes = json.load(open(recipes_path))
        self.hof_fitness_sum = False

        if githash is not None:
            if self.etype is None:
                meta = runlog[str(githash)][str(seed)]
                self.etype = meta['etype']
                self.notes = meta['notes']
        elif self.etypetest is not None:
            self.etype = etypetest
        else:
            print('No test etype or githash defined, nothing to do!')
            exit(1)

        # load from checkpoint, get parameters
        if githash is not None:

            self.checkpoints_dir = 'checkpoints/run.%s' % githash
            cp_filename = os.path.join(
                self.checkpoints_dir, 'checkpoint_%d.pkl' % seed)

            if os.path.isfile(cp_filename):
                pass
            else:
                print('No checkpoint file available, run optimization first')
                exit(1)

            cp = pickle.load(open(cp_filename, "r"))
            hof = cp['halloffame']

            if self.has_etest:
                evaluator = setup.evaluator.create(etype=self.etype,
                                        altmorph=self.altmorph,
                                        etypetest=self.etypetest)
            else:
                evaluator = setup.evaluator.create(etype=self.etype,
                                        altmorph=self.altmorph)
            self.parameters = evaluator.param_dict(hof[rank])
            self.cp = cp
            self.hof_fitness_sum = sum(hof[rank].fitness.values)
            #self.hof_objectives = evaluator.objective_dict(hof[rank].fitness.values)


    def get_name(self):

        altmorph = self.altmorph
        if altmorph is not None:
            if isinstance(altmorph, (list)):
                altmorph = altmorph[0][1]
            altmorph = str(os.path.basename(altmorph)).split(".")[0]

        # define output label
        if self.etypetest is not None:
            reportname = 'etype:%s-test ' % self.etypetest
            label = reportname
            report_elem = OrderedDict()
            if altmorph is not None:
                report_elem['altmorph'] = altmorph
        else:
            reportname = ''
            label = ''
            if altmorph is None:
                report_elem = OrderedDict([('etype',self.etype), ('githash',self.githash)
                                        ,('seed',self.seed), ('rank',self.rank)])
            else:
                report_elem = OrderedDict([('etype',self.etype), ('githash',self.githash)
                                        ,('altmorph',altmorph)
                                        ,('seed',self.seed), ('rank',self.rank)])


        # stimulation options
        if self.usethreshold:
            report_elem['usethreshold'] = 'True'

        if self.stochdet is not None:
            report_elem['stochdet'] = str(self.stochdet)


        for elem_name, elem_val in report_elem.iteritems():
            if elem_name in self.grouping:
                reportname += '%s:%s ' % (elem_name,elem_val)
            else:
                label += '%s:%s ' % (elem_name,elem_val)

        return reportname.rstrip(), label.rstrip()


    def plot_evolution(self, figs, color='b'):

        if hasattr(self, 'cp'):
            reportname, label = self.get_name()
            evol_fig = plot.evolution(self.cp['logbook'], figs=figs,
                                color=color, reportname=reportname)
            plt.show(block=False)


    def plot_diversity(self, figs, color='b'):

        if hasattr(self, 'cp'):
            reportname, label = self.get_name()
            if self.has_etest:
                evaluator = self.setup.evaluator.create(etype=self.etype,
                                    altmorph=self.altmorph,
                                    etypetest=self.etypetest)
            else:
                evaluator = self.setup.evaluator.create(etype=self.etype,
                                    altmorph=self.altmorph)

            evol_fig = plot.diversity(self.cp, evaluator=evaluator,
                                color=color, figs=figs,
                                reportname=reportname)
            plt.show(block=False)


    def do_model_export(self):

        if hasattr(self, 'cp') and (self.altmorph is None):
            # generate hoc

            if self.has_etest:
                evaluator = self.setup.evaluator.create(etype=self.etype,
                                    etypetest=self.etypetest)
            else:
                evaluator = self.setup.evaluator.create(etype=self.etype)

            # also make it compatible to old single-morph evaluators
            if hasattr(evaluator, 'evaluators'):
                evaluators = evaluator.evaluators
            else:
                evaluators = [evaluator]

            for i, evl in enumerate(evaluators):
                hoccode = evl.cell_model.create_hoc(param_values=self.parameters)
                hoc_path = os.path.join(self.checkpoints_dir, "%s_%s_%s_%s.hoc" % (i, self.etype, self.githash, self.seed))
                with open(hoc_path, "w") as f:
                    f.write(hoccode)

            #generating parameter definition for this individual
            params_path = self.recipes[self.etype]['params']
            with open(params_path) as params_file:
                definitions = json.load(
                    params_file,
                    object_pairs_hook=OrderedDict)
            params_definitions = definitions["parameters"]

            for param_name, param_value in self.parameters.iteritems():
                name = param_name.split(".")[0]
                location = param_name.split(".")[1]
                for param in params_definitions[location]:
                    if name == param["name"]:
                        param["val"] = param_value
                        if "test" in param:
                            del param["test"]

            path = os.path.join(self.checkpoints_dir, "%s_%s_%s.json" % (self.etype, self.githash, self.seed))
            s = json.dumps(definitions, indent=2)
            s = self.collapse_json(s, indent=6)
            with open(path, "w") as f:
                f.write(s)


    def create_cell_model(self, evaluator):

        self.sim = evaluator.sim

        # only use alternative mode, e.g. hoc if etypetest is given
        if self.etypetest is not None:

            # generate parameters
            parameters = {}

            params_path = self.recipes[self.etypetest]['params']
            with open(params_path) as params_file:
                definitions = json.load(
                    params_file,
                    object_pairs_hook=OrderedDict)
            params_definitions = definitions["parameters"]

            if "__comment" in params_definitions:
                del params_definitions["__comment"]

            for sectionlist, params in params_definitions.iteritems():
                for param_config in params:
                    param_name = param_config["name"]
                    if isinstance(param_config["val"], (list, tuple)):
                        test = param_config["test"]
                        parameters['%s.%s' % (param_name, sectionlist)] = test

            self.parameters = parameters

            for evl in self.evaluators:

                morph = evl.cell_model.morphology.morphology_path
                if self.hoc:
                    evl.sim.neuron.h.celsius = 34
                    evl.sim.neuron.h.v_init = -80
                    evl.cell_model = ephys.models.HocCellModel(
                                        'hoc', morph, self.hoc)

                elif self.oldhoc:
                    import setup.hocmodel
                    evl.sim.neuron.h.celsius = 34
                    evl.sim.neuron.h.v_init = -80
                    evl.cell_model = setup.hocmodel.HocOldModel(
                                        "oldhocmodel", "./morphologies", self.oldhoc)


    def sim_objectives(self, stochdet=None, usethreshold=False):

        self.stochdet = stochdet
        self.usethreshold = usethreshold

        if self.has_etest:
            evaluator = self.setup.evaluator.create(etype=self.etype,
                                stochkv_det=stochdet,
                                usethreshold=self.usethreshold,
                                altmorph=self.altmorph,
                                etypetest=self.etypetest)
        else:
            evaluator = self.setup.evaluator.create(etype=self.etype,
                                stochkv_det=stochdet,
                                usethreshold=self.usethreshold,
                                altmorph=self.altmorph)



        # also make it compatible to old single-morph evaluators
        if hasattr(evaluator, 'evaluators'):
            self.evaluators = evaluator.evaluators
        else:
            self.evaluators = [evaluator]

        self.create_cell_model(evaluator)


        start_time_all = time.time()

        #features = []
        features = {}
        responses = {}
        objectives = {}

        for evl in self.evaluators:

            fitness_protocols = evl.fitness_protocols

            for protocol in fitness_protocols.values():
                start_time = time.time()
                response = evl.run_protocol(protocol,
                    param_values=self.parameters)
                responses.update(response)
                logger.info(" Ran protocol in %f seconds",
                                time.time() - start_time)

            for obj in evl.fitness_calculator.objectives:
                for feature in obj.features:
                    if feature.calculate_feature(responses) is not None:
                        features[feature.name] = round(feature.calculate_feature(responses),2)
                    else:
                        features[feature.name] = feature.calculate_feature(responses)

                        #features.append(feature.calculate_feature(responses))

            score = evl.fitness_calculator.calculate_scores(responses)
            objectives.update(score)

        self.responses = responses
        self.objectives = objectives
        self.features = features

        logger.info(" Full evaluation took %f seconds",
                        time.time() - start_time_all)


    def plot_obj(self, figs, color='b', split_sim=1):

        reportname, label = self.get_name()
        responses_fig = plot.responses(self.responses, figs=figs,
                        color=color, cols=split_sim,
                        reportname=reportname,
                        label=label)
        plt.show(block=False)

        objectives_fig = plot.objectives(self.objectives, self.features, figs=figs,
                            color=color, reportname=reportname,
                            label=label)

        if self.hof_fitness_sum:
            objsum = sum(self.objectives.values())
            msg = ("Sum of scores from hof: %s, from simulation: %s" %
                            (self.hof_fitness_sum, objsum) )
            logger.info(" " + msg)

            if (self.usethreshold is False) and (self.stochdet is None):
                # check if equal to optimization
                if abs(objsum-self.hof_fitness_sum) > 0.5 * 10**(-6):
                    color = 'red'
                    prefix = 'ERROR:'
                else:
                    color = 'green'
                    prefix = 'GOOD:'

                msg = (prefix + " Sum of scores " +
                        "from HOF: %0.2f, from simulation: %0.2f diff: %0.4f" %
                        (self.hof_fitness_sum, objsum, abs(objsum-self.hof_fitness_sum)) )

                objectives_fig.text(0.05, 0.5, msg,
                         fontsize=15, color=color,
                         ha='center', va='center',
                         rotation=90)

            plt.show(block=False)

        return reportname


    def sim_plot_obj(self, stochdet=None, usethreshold=False,
                    figs={}, color='b', split_sim=1):

        self.sim_objectives(stochdet=stochdet, usethreshold=usethreshold)
        self.plot_obj(figs=figs, color=color, split_sim=split_sim)


    def save_pdf(self, reportname, figs, subdir='figures'):
        reportfile = reportname.rstrip().replace (" ", "_").replace (":", "_")
        report_path = os.path.join(self.figpath, reportfile + '.pdf')
        pdf_pages = PdfPages(report_path)
        for figname, fig in figs.iteritems():
            pdf_pages.savefig(fig['fig'])
        pdf_pages.close()


    def add_to_final(self, branch, emodel):

        #path_rundir = os.path.dirname(os.path.realpath(__file__))
        #path_rundir = os.path.join(path_rundir, 'run', githash)
        
        with open(self.path_final) as f:
            final = json.load(f,
                object_pairs_hook=OrderedDict)
        print self.githash
        print branch
        # delete old entry
        final[emodel] = OrderedDict()

        final[emodel]["branch"] = str(branch)
        final[emodel]["githash"] = self.githash
        final[emodel]["seed"] = self.seed
        final[emodel]["rank"] = self.rank
        final[emodel]["notes"] = self.notes
        final[emodel]["params"] = self.parameters
        final[emodel]["fitness"] = self.objectives
        final[emodel]["score"] = sum(self.objectives.values())
        final[emodel]["morph_path"] = os.path.relpath(self.evaluators[0].cell_model.morphology.morphology_path)
        print final[emodel]
        s = json.dumps(final, indent=2)
        s = self.collapse_json(s, indent=4)
        with open(self.path_final, "w") as f:
            f.write(s)
