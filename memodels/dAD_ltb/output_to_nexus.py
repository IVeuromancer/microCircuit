import os
import glob
import requests
from entity_management.simulation import (EModelScript, EModelRelease,
                                          ModelReleaseIndex, Distribution, BluePyOptRun)
from input_from_nexus import patch_run, run_name, use_nexus_for_thalamus
import logging
logger = logging.getLogger()

def check_token_validity(token):
    url = 'https://bbp-nexus.epfl.ch/staging/v0/oauth2/userinfo'
    r = requests.get(url, headers={'authorization': token})
    if token is None:
        raise Exception('Environment variable NEXUS_TOKEN is empty. It should contain a Nexus token. You can get one by going to https://bbp-nexus.epfl.ch/staging/explorer/ and clicking the "Copy token" button')
    if r.status_code == 200:
        return
    print("r.text: {}".format(r.text))
    if r.status_code == 500:
        raise Exception('{} is returning an Error 500. Nexus is probably down. '
                        'Try again later')
    if r.status_code == 401:
        raise Exception('Error 401: your token has expired or you are not authorized.\n'
                        'Suggestion: try renewing the token stored in the environment variable: NEXUS_TOKEN')
    raise Exception('Received error code: {}, {}'.format(r.status_code, r.text))

check_token_validity(os.getenv('NEXUS_TOKEN'))


def upload_hoc_files_in_dir(folder, token):
    for filepath in glob.glob(folder+'/*.hoc'):
        upload_hoc_file(filepath, token)

def upload_emodel_release(filepath, githash, seed, token):
    if not use_nexus_for_thalamus() or githash is None:
        return

    name = 'thalamus BluePyOpt emodel release {} #{}'.format(githash, seed)
    logger.warning('Uploading: {}'.format(name))

    def create_release():
        path = '{}/thalamus/run/{}'.format(filepath, githash)
        index_name = 'BluePyOpt index {}, seed: {}'.format(githash, seed)
        index = ModelReleaseIndex.find_unique(
            name=index_name,
            on_no_result=lambda: ModelReleaseIndex(
                name=index_name,
                distribution=[Distribution(
                    downloadURL='file://{}/final.json'.format(path),
                    mediaType='application/json')]
            ).publish()
        )

        model = EModelRelease(
            name=name,
            distribution=[
                Distribution(downloadURL='file://{}/checkpoints/run.{}'.format(path, githash))
            ],
            emodelIndex=index
        )

        return model.publish(use_auth=token)

    release = EModelRelease.find_unique(name=name, on_no_result=create_release)
    run = BluePyOptRun.find_unique(name=run_name(githash, seed), throw=True)

    patch_run(run.id,
              'hasOutput',
              release)




def upload_hoc_file(filepath, token):
    name = os.path.basename(filepath)
    emodel_script = EModelScript(name=name).publish(use_auth=token)
    with open(filepath) as f:
        emodel_script.attach(name, f, 'application/neuron-hoc', use_auth=token)



if __name__=='__main__':
    upload_emodel_release('/home/bcoste/workspace/nexus/optimization-mouse-nexus',
                          '043a017',
                          os.getenv('NEXUS_TOKEN'))
