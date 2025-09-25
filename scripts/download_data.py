"""
A script to download the from Thanh et al. (2017).

While we have made an effort to make sure that the data retrieved is correct,
we do not own and do not have the rights to it, and do not guarantee its correctness.
"""
from pathlib import Path
from typing import Any
import urllib.request
import urllib.error
import urllib.parse
import hashlib
import zipfile


BASE_PATH = Path(__file__).parent.parent

MODELS_URL = r'https://www.cosbi.eu/wp-content/uploads/2021/11/Collection-of-models-RSSA.zip'
MODELS_ZIP_PATH = BASE_PATH / 'data/models.zip'
MODELS_DIR = BASE_PATH / 'data/models'
MODEL_PATHS = [MODELS_DIR / 'B cell antigen receptor signaling', MODELS_DIR / 'FceRI']


MODELS_HASH = b'\xad(\xc6\x17\xd9\x9b\xc4\x9e\x80\xdd\xcd\x04M~\x85\x08S\xcd>\xe5s\xceF\xc4"|\xa4\xf1\x96\xd2(\xb1\x7f\x04-?Z\xb1C\xf12p\x82\x9c\xba\n\xbc\xab\xe9\x98\xe3ktv\x9e\x9e\x13D\xb1v\xb3\x1b\xbd\x12'
REPO_URL = r"https://raw.githubusercontent.com/debraj86/Algorithms/226275e08e8931fe0236b9832d280b51cdbc9daa"
NETWORK_NAMES = ['B cell receptor signaling network', 'Fceri signaling network']

POP_NAME = 'conc.txt'
REACTION_RATE_NAME = 'k.txt'
STOICHIOMETRY_MATRIX_NAME = 'xyz.txt'

def file_hash(path: Path, hash_alg:Any=hashlib.sha512) -> bytes:
    return hash_alg(path.open('rb').read()).digest()


def get_thanh_models():
    """
    Gets the model data from Thanh et al. (2017).
    """
    (BASE_PATH / 'data').mkdir(exist_ok=True)

    if all(path.exists() for path in MODEL_PATHS):
        print('Models directory found.')
        return
    
    if not MODELS_ZIP_PATH.exists() or file_hash(MODELS_ZIP_PATH) != MODELS_HASH:
        print('Downloading the models.')
        headers = {'User-Agent': 'Mozilla/5.0'}
        request = urllib.request.Request(MODELS_URL, headers=headers)
        try:
            with urllib.request.urlopen(request) as response:
                # Read the content directly and write it to a local file.
                with open(MODELS_ZIP_PATH, 'wb') as out_file:
                    _ = out_file.write(response.read())
        except urllib.error.HTTPError as e:
            raise RuntimeError(f'Failed to download the models from {MODELS_URL} with response: {e.code} - {e.reason}. '
            'If the problem persists, please download the model zip manually and place it in the data directory.')

    if not MODELS_ZIP_PATH.exists():
        raise RuntimeError(f'Failed to download the models from {MODELS_URL}. Please download the zip file and place it in the data directory.')

    if file_hash(MODELS_ZIP_PATH) != MODELS_HASH:
        raise RuntimeError(f'The model zip downloaded from {MODELS_URL} is not the expected file. '
        'If you believe it stores the correct files, please extract them into data/model.')

    print('Extracting the models.')
    MODELS_DIR.mkdir(exist_ok=True)
    zipfile.ZipFile(MODELS_ZIP_PATH).extractall(BASE_PATH / 'data')

def get_file(relative_path: str, base_path:str=REPO_URL) -> str:
    file_url = f'{base_path}/{urllib.parse.quote(relative_path)}'
    with urllib.request.urlopen(file_url) as response:
        # Read the content as bytes and decode it into a string.
        return response.read().decode('utf-8')

def get_blsssa_network(network: str, out_file: Path):
    """
    Downloads the BlSSSA network, and converts it to the Tau-Splitting format.

    The BlSSSA networks do not allow shared inputs and outputs, 
    cancelling them out and removing them from the stoichiometry.
    """
    initial_state = [int(row) for row in get_file(f'{network}/{POP_NAME}').split()]
    reaction_rates = [float(row) for row in get_file(f'{network}/{REACTION_RATE_NAME}').split()]
    reactions = [[int(x.strip()) for x in row.split(',') if x.strip()] for row in get_file(f'{network}/{STOICHIOMETRY_MATRIX_NAME}').split('\n')]

    out_file.parent.mkdir(exist_ok=True)

    with open(out_file, 'w') as f:
        for idx, initial_count in enumerate(initial_state):
            f.write(f'S{idx + 1} = {initial_count}\n')
        

        for r_idx, rate in enumerate(reaction_rates):
            inputs = []
            outputs = []
            for c_idx, row in enumerate(reactions):
                if not row:
                    continue
                if row[r_idx] < 0:
                    if row[r_idx] < -1:
                        inputs.append(f'{-row[r_idx]}S{c_idx + 1}')
                    else:
                        inputs.append(f'S{c_idx + 1}')
                if row[r_idx] > 0:
                    if row[r_idx] > 1:
                        outputs.append(f'{row[r_idx]}S{c_idx + 1}')
                    else:
                        outputs.append(f'S{c_idx + 1}')
            
            f.write(f'{' + '.join(inputs)} -> {' + '.join(outputs)}, {rate}\n')
    

if __name__ == '__main__':
    get_thanh_models()
    get_blsssa_network('B cell receptor signaling network', MODELS_DIR / 'blsssa_bcr.txt')
    get_blsssa_network('Fceri signaling network', MODELS_DIR / 'blsssa_fceri.txt')


