from dataclasses import dataclass
from pathlib import Path


BASE_PATH = Path(__file__).parent.parent
MODELS_DIR = BASE_PATH / 'data/models'

MODEL_PATHS = [MODELS_DIR / 'B cell antigen receptor signaling', MODELS_DIR / 'FceRI']

INITIAL_STATE_FILENAME = 'conc.txt'
RATE_FILENAME = 'k.txt'
STOI_FILENAME = 'xyz.txt'
INP_PROD_FILENAME = 'reaction_inputs.txt'

@dataclass
class Reaction:
    inputs: dict[str, int]
    outputs: dict[str, int]
    rate: float

def parse_network(initial_state_file: Path, reactions_file: Path) -> tuple[dict[str, int], list[Reaction]]:
    initial_state = {}
    reactions = []
    for line in initial_state_file.open().readlines():  
        parts = line.strip().split('=')
        name = parts[0].strip()
        amount = int(parts[1].strip())
        initial_state[name] = amount

    for line in reactions_file.open().readlines():
        line, rate = line.split(',')

        input_count = {}
        output_count = {}

        rate = float(rate.strip())
        inputs, outputs = line.split('->')

        for inp in  inputs.split('+'):
            input_count[inp.strip()] = input_count.get(inp.strip(), 0) + 1
        for outp in  outputs.split('+'):
            output_count[outp.strip()] = output_count.get(outp.strip(), 0) + 1
        
        reactions.append(Reaction(input_count, output_count, rate))

    return initial_state, reactions

def print_network(initial_state: dict[str, int], reactions: list[Reaction], out_dir: Path):
    """
    Prints the chemical reaction network data in the BlSSSA format.
    """
    if not out_dir.exists():
        out_dir.mkdir()

    reactant_idx = {name: idx for idx, name in enumerate(initial_state)}

    initial_state_vec = [0] * len(initial_state)
    for name, idx in reactant_idx.items():
        initial_state_vec[idx] = initial_state[name]

    with open(out_dir / INITIAL_STATE_FILENAME, 'w') as f:
        f.write('\n'.join(map(str, initial_state_vec)))
    
    with open(out_dir / RATE_FILENAME, 'w') as f:
        f.write('\n'.join(str(r.rate) for r in reactions))

    stoi_matrix = [[0]*len(reactions) for _ in initial_state]
    input_matrix = [[0]*len(reactions) for _ in initial_state]

    for r_idx, reaction in enumerate(reactions):
        for inp, count in reaction.inputs.items():
            stoi_matrix[reactant_idx[inp]][r_idx] -= count
            input_matrix[reactant_idx[inp]][r_idx] += count
        for outp, count in reaction.outputs.items():
            stoi_matrix[reactant_idx[outp]][r_idx] += count

    with open(out_dir / STOI_FILENAME, 'w') as f:
        for row in stoi_matrix:
            f.write(f'{','.join(map(str, row))},\n')

    with open(out_dir / INP_PROD_FILENAME, 'w') as f:
        for row in input_matrix:
            f.write(f'{','.join(map(str, row))},\n')



print_network(*parse_network(MODEL_PATHS[0] / 'BCR_pop.txt', MODEL_PATHS[0] / 'BCR_rxn.txt'), MODEL_PATHS[0] / 'blsssa')
print_network(*parse_network(MODEL_PATHS[1] / 'BCR_pop.txt', MODEL_PATHS[1] / 'BCR_rxn.txt'), MODEL_PATHS[1] / 'blsssa')
