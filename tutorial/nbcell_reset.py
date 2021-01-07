import json
from pathlib import Path


def reset(nb_file):
    with open(nb_file, 'r') as f_in:
        doc = json.load(f_in)

    count = 1

    for cell in doc['cells']:
        if 'execution_count' not in cell:
            continue

        cell['execution_count'] = count

        for o in cell.get('outputs', []):
            if 'execution_count' in o:
                o['execution_count'] = count

        count += 1

    with open(nb_file, 'w') as f_out:
        json.dump(doc, f_out, indent=1)


if __name__ == "__main__":
    for nb in Path().cwd().glob("*.ipynb"):
        reset(nb)
