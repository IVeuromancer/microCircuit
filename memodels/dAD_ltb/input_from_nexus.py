import json
import os
import requests
import tempfile
from time import sleep

import logging
logger = logging.getLogger(__name__)


from entity_management import from_url
from entity_management.morphology import ReconstructedWholeBrainCell, ReconstructedPatchedCell
from entity_management.simulation import BluePyOptRun, Configuration, SubCellularModelScript, Distribution, IonChannelMechanismRelease, BluePyEfeFeatures

run_url = 'https://bbp-nexus.epfl.ch/staging/v0/data/thalamusproject/simulation/bluepyoptrun/v0.1.8'
def use_nexus_for_thalamus():
    return os.getenv('USE_NEXUS_FOR_THALAMUS', False)


def recipes_from_nexus():
    if not use_nexus_for_thalamus():
         return './config/recipes/recipes.json'
    return Configuration.find_unique(name="Thalamus recipe", throw=True).distribution[0].downloadURL[7:]

def _get_configuration(etype, stuff, filter_on_name):
    with open(recipes_from_nexus()) as f:
        recipe = json.load(f)

    if not use_nexus_for_thalamus():
        return None, os.path.join('.',
                            recipe[etype][stuff])

    name = filter_on_name(os.path.basename(recipe[etype][stuff]))
    stuff = Configuration.find_unique(name=name)
    if not stuff:
        raise Exception('Configuration: ({}) not found'.format(name))

    return stuff, stuff.distribution[0].downloadURL[7:]


def patch_run(run_url, key, value):
    if not use_nexus_for_thalamus() or run_url is None:
        return
    logging.warning('Adding {} to BluePyOptRun'.format(key))
    return  from_url(run_url).evolve(**{key: value}).publish()


_CURRENT_RUN=None

def run_name(githash, seed):
    return "BluePyOptRun: {}-{}".format(githash, seed)

def register_current_run(githash, seed):
    if not use_nexus_for_thalamus() or githash is None:
        return
    global _CURRENT_RUN
    name = run_name(githash, seed)
    logging.info("Registering run: {}".format(name))
    _CURRENT_RUN = BluePyOptRun.find_unique(
        name=name,
        on_no_result=lambda: BluePyOptRun(
            name=name,
            gitHash=githash,
            inputMechanisms=IonChannelMechanismRelease.find_unique(name='Thalamus release')
        ).publish(),
        poll_until_exists=True
    ).id

    print("_CURRENT_RUN: {}".format(_CURRENT_RUN))

    for _ in range(12):
        run = BluePyOptRun.find_unique(name=name)
        if run:
            break
        print('Waiting for the run to be registered...')
        sleep(10)
    if not run:
        raise Exception('Timeout, cannot retrieve run: {}. '.format(name))



def protocol_from_nexus(etype):
    protocol, path = _get_configuration(etype, 'protocol', lambda name: "Thalamus protocol: {}".format(name))
    patch_run(_CURRENT_RUN, "bluePyOptProtocol", protocol)
    return path

def params_from_nexus(etype):
    params, path = _get_configuration(etype, 'params', lambda name: "Thalamus params: {}".format(name))
    patch_run(_CURRENT_RUN, "bluePyOptParameters", params)
    return path

def morphology_from_nexus(etype, altmorph):
    if altmorph is None:
        with open(recipes_from_nexus()) as f:
            recipe = json.load(f)
        morph_path = os.path.join(os.path.join(recipe[etype]['morph_path'], recipe[etype]['morphology']))
    else:
        morph_path = altmorph

    if use_nexus_for_thalamus():
        name = recipe[etype]['morphology'][:-4]
        morpho = (ReconstructedWholeBrainCell.find_unique(name=name) or
                  ReconstructedPatchedCell.find_unique(name=name))
        patch_run(_CURRENT_RUN, "morphology", morpho)

    return morph_path


def features_from_nexus(etype):
    with open(recipes_from_nexus()) as f:
        recipe = json.load(f)
    if not use_nexus_for_thalamus():
        return os.path.join('.',
                            recipe[etype]['features'])


    print("recipe[etype]['features']: {}".format(recipe[etype]['features']))
    feature = BluePyEfeFeatures.find_unique(name=recipe[etype]['features'])
    run = patch_run(_CURRENT_RUN, "experimentalFeatures", feature)
    filename = tempfile.NamedTemporaryFile().name
    feature.features.download(os.path.join('/tmp', filename))
    path = os.path.join('/tmp', filename)

    return path
