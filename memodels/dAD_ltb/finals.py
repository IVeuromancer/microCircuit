#!/usr/bin/env python

import sys
import traceback
import os
import sh
import json
import matplotlib
from output_to_nexus import upload_emodel_release
import collections
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

import zlib
import zipfile
import gzip
from setup.tools import collapse_json
import subprocess


def str2bool(v):
    if str(v).lower() in ("yes", "true", "t", "1"):
        return True
    elif str(v).lower() in ("no", "false", "f", "0"):
        return False
    elif str(v).lower() in ("none"):
        return None
    else:
        return v

import argparse
parser = argparse.ArgumentParser(description='update')
parser.add_argument('--writerundir', action="store_true")
parser.add_argument('--addfiles', action="store_true")
parser.add_argument('--report', action="store_true")
parser.add_argument('--patch', action="store_true")
parser.add_argument('--run', action="store_true")
parser.add_argument('--branch', type=str, required=False, default="None")
parser.add_argument('--release', type=str, required=False, default="False")

args, unknown = parser.parse_known_args()

writerundir = args.writerundir
addfiles = args.addfiles
doreport = args.report
patch = args.patch
run = args.run
argbranch = args.branch
release = str2bool(args.release)

main_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(main_path)

if os.getenv('USEIPYP') == '1':
    logfile = os.path.join(main_path, 'logs/report.log')
    logging.basicConfig(filename=logfile)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.info(' Using ipyparallel')

    matplotlib.use('pdf')
    from ipyparallel import Client
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
    view = rc.load_balanced_view()

    map_function = view.map_sync
    #map_function = view.imap

    do_plot = False

else:
    logging.basicConfig(stream=sys.stdout)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.info(' Using map')

    map_function = map
    matplotlib.use('TkAgg')


import report


def make_release(releasehash=None):
    """Main"""

    if releasehash is None:
        print("Making new release commit")
        print sh.git.commit('--allow-empty', '-a', m='Latest release commit!')
        releasehash = sh.git('rev-parse', '--short', 'HEAD').rstrip()

    subprocess.call("git archive --format=tar --prefix=%s/ %s emodel_etype_map.json | (cd ./releases/ && tar xf -)" % (releasehash, releasehash), shell=True)
    subprocess.call("git archive --format=tar --prefix=%s/ %s final.json | (cd ./releases/ && tar xf -)" % (releasehash, releasehash), shell=True)
    featureextractpath = "/gpfs/bbp.cscs.ch/project/proj38/singlecell/features"

    emap = json.load(open("./releases/%s/emodel_etype_map.json" % releasehash, 'r'),
        object_pairs_hook=collections.OrderedDict)

    final = json.load(open("./releases/%s/final.json" % releasehash, 'r'),
        object_pairs_hook=collections.OrderedDict)

    output = collections.OrderedDict()
    output_scores = collections.OrderedDict()

    for emname, emodel in emap.iteritems():

        if "mtype" in emap:
            mtype = emodel["mtype"]
        else:
            mtype = ".*"

        etype = emodel["etype"]

        if etype not in output:
            output[etype] = []#collections.OrderedDict()
            output_scores[etype] = []

        #if mtype not in output[etype]:
        #    output[etype][mtype] = []

        d = collections.OrderedDict()
        s = collections.OrderedDict()

        layer = emodel["layer"]

        print etype, mtype, layer, emname

        d["mtype"] = mtype
        d["layer"] = layer
        d["emname"] = emname
        d["features"] = collections.OrderedDict()

        s["mtype"] = mtype
        s["layer"] = layer
        s["emname"] = emname
        s["features"] = collections.OrderedDict()

        mfinal = final[emname]
        branch = mfinal['branch']
        fitness = mfinal['fitness']

        if True:
            recipepath = 'config/recipes/recipes.json'
            sh.git.archive('--output=recipes.json.zip', branch, recipepath)
            with zipfile.ZipFile("recipes.json.zip", 'r') as z:
                content = z.read(recipepath)
                recipes = json.loads(content, object_pairs_hook=collections.OrderedDict)

            try:
                os.remove('recipes.json.zip')
            except OSError:
                pass

            featurepath = recipes[emname]['features']
            sh.git.archive('--output=features.json.zip', branch, featurepath)
            with zipfile.ZipFile("features.json.zip", 'r') as z:
                content = z.read(featurepath)
                features = json.loads(content, object_pairs_hook=collections.OrderedDict)
            os.remove('features.json.zip')

            try:
                os.remove('features.json.zip')
            except OSError:
                pass

            githash = None
            source = None

            for expname, experiment in features.iteritems():
                d["features"][expname] = collections.OrderedDict()
                s["features"][expname] = collections.OrderedDict()

                for unitname, unit in experiment.iteritems():
                    d["features"][expname][unitname] = collections.OrderedDict()
                    s["features"][expname][unitname] = collections.OrderedDict()

                    for feature in unit:

                        # collect scores, save them

                        if True:
                            featname = feature["feature"]

                            print expname, unitname, featname

                            d["features"][expname][unitname][featname] = collections.OrderedDict()
                            feat = d["features"][expname][unitname][featname]

                            s["features"][expname][unitname][featname] = collections.OrderedDict()
                            scores = s["features"][expname][unitname][featname]

                            # report on scores
                            prefix = fitness.keys()[0].split('.')[0]
                            featurestring = prefix + "." + expname + "." + unitname + "." + featname
                            fit = fitness[featurestring]

                            # remove weight used for fitting!
                            if "weight" in feature:
                                fit = fit / feature["weight"]

                            scores['score'] = fit

                            if ("__comment" in feature) and isinstance(feature["__comment"], basestring):
                                scores['comment'] = feature["__comment"]
                                feat['comment'] = feature["__comment"]

                            elif (featname == "bpo_holding_current") or (featname == "bpo_threshold_current"):
                                scores['comment'] = "As set by experimenter"
                                feat['comment'] = "As set by experimenter"

                            elif "fid" not in feature:
                                scores['comment'] = "Anecdotal value"
                                feat['comment'] = "Anecdotal value"

                            mean = feature["val"][0]
                            std = feature["val"][1]

                            feat['mean'] = mean
                            feat['std'] = std

                            # only add if recordings were made!
                            if "fid" in feature:

                                fid = feature["fid"]

                                if ( (githash != feature["__comment"]["meta"]["version"])
                                    and
                                    (source != feature["__comment"]["source"]) ):

                                    githash = feature["__comment"]["meta"]["version"]
                                    source = feature["__comment"]["source"]

                                    featurespath = source + '/features_sources.json.gz'
                                    sh.git.archive(
                                        '--output=features_sources.json.gz.zip',
                                        '--remote=' + featureextractpath,
                                        githash,
                                        featurespath)

                                    with zipfile.ZipFile("features_sources.json.gz.zip", 'r') as z:
                                        content = z.read(featurespath)

                                        features_sources = json.loads(
                                                zlib.decompress(content, 16+zlib.MAX_WBITS),
                                                object_pairs_hook=collections.OrderedDict)

                                    try:
                                        os.remove('features_sources.json.gz.zip')
                                    except OSError:
                                        pass

                                    allfeat = collections.OrderedDict(
                                        (f["fid"], f)
                                        for en, exp in features_sources.items()
                                        for un, unit in exp.items()
                                        for f in unit
                                    )


                                raw = allfeat[fid]['raw']

                                for l1 in raw:
                                    for l2 in l1:
                                        if (l2['t_unit'] == ""):
                                            l2['t_unit'] = "s"

                                feat['raw'] = raw



        #                 except:
        #                     pass
        # except:
        #     pass

        output[etype].append(d)
        output_scores[etype].append(s)

    s = json.dumps(output_scores, indent=2)
    s = collapse_json(s, indent=14)
    with open("releases/%s/scores_collection.json" % releasehash, "w") as f:
        f.write(s)

    s = json.dumps(output, indent=2)
    s = collapse_json(s, indent=14)
    with gzip.open("releases/%s/features_collection.json.gz" % releasehash, "wb") as f:
        f.write(s)

    print("Written release %s", releasehash)

    return releasehash

def main():
    """Main"""

    if release is not False:

        releasehash = make_release(releasehash=release)

        recipepath = 'final.json'
        sh.git.archive('--output=final.json.zip', releasehash, recipepath)
        with zipfile.ZipFile("final.json.zip", 'r') as z:
            content = z.read(recipepath)
            finals = json.loads(content, object_pairs_hook=collections.OrderedDict)
        os.remove('final.json.zip')
        figpath = 'releases/%s/report' % releasehash

        print("Release hash given!")

    else:
        finals = json.load(open('final.json'))
        figpath = 'figures'

    for name, final in finals.iteritems():
        print('calling upload for: {}'.format(name))
        upload_emodel_release(main_path,
                              final["githash"],
                              final["seed"],
                              os.getenv('NEXUS_TOKEN'))

    if doreport:

        args_list = []
        for name, final in finals.iteritems():
            githash = final["githash"]
            seed = final["seed"]
            rundir = "run/%s" % githash

            args = {'githash':githash,
                    'etype':name,
                    'seed':seed,
                    'main_path':main_path,
                    'rundir':rundir,
                    'figpath':figpath
                    }
            print "args:", args
            args_list.append(args)

        map_function(report.report_isolated, args_list)

    elif patch:

        message = "patching setup/evaluator.py"

        #print sh.git.add('-A')
        print sh.git.commit('--allow-empty', '-a', m='committing before patching')

        for name, final in finals.iteritems():

            print name

            sh.git.checkout(final["branch"])

            os.system('git diff HEAD..master setup/evaluator.py > file.patch')
            #patch = sh.git.diff('HEAD..master', 'setup/evaluator.py') # , _out="files.list"

            os.system('filterdiff --lines=480- file.patch > file_.patch')
            os.remove('file.patch')

            patchstr = open('file_.patch', 'r').read()

            print patchstr

            if len(patchstr) > 0:
                if run:
                    print sh.git.apply('file_.patch')
                    os.remove('file_.patch')
                    print sh.git.commit('-a', m=message)
                else:
                    print sh.git.apply('--check', 'file_.patch')
                    os.remove('file_.patch')
            else:
                os.remove('file_.patch')

        sh.git.checkout('master')


    elif writerundir:

        for name, final in finals.iteritems():
           branch = final["branch"]
           githash = final["githash"]
           command = "git archive --format=tar --prefix=%s/ %s | (cd ./run/ && tar xf -)" % (githash, branch)
           result = os.system(command)
           print "Updating", branch, githash, result

    elif addfiles:
        message = "updating all finals to inlcude base seed based on params, bugfixes"
        files = [
                    # "./config/protocols/bAC.json",
                    # "./config/protocols/bNAC.json",
                    # "./config/protocols/cACint.json",
                    # "./config/protocols/cNAC.json",
                    # "./config/protocols/dNAC.json",
                    # "./config/protocols/cADpyr_L4PC.json",
                    # "./config/protocols/cADpyr_L5PC.json",
                    # "./config/protocols/cADpyr_L6UTPC.json",
                    # "./config/protocols/cADpyr_L23PC.json"
                    # "./config/recipes/recipes.json"
                    #"./setup/template.py",
                    "./setup/evaluator.py",
                    #"./setup/protocols.py",
                    #"./opt_model.py"
                ]

        print sh.git.add('-A')
        print sh.git.commit('-a', m='committing before running replacing files')

        if (argbranch is not "None"):
            sh.git.checkout(argbranch)

            for f in files:
                print sh.git.checkout('master', f)

            print sh.git.commit(m=message)
            print sh.git.status()

        else:
            for name, final in finals.iteritems():
                sh.git.checkout(final["branch"])

                for f in files:
                    print sh.git.checkout('master', f)

                print sh.git.commit(m=message)
                print sh.git.status()

        sh.git.checkout('master')

if __name__ == '__main__':
    try:
        main()
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
