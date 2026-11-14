"""Frozen E3 design and stable stream identities, independent across configs."""
from pathlib import Path
import hashlib
import itertools
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]


def read_config(path=None):
    config = yaml.safe_load(Path(path or ROOT/'experiments/e3/config.yaml').read_text())
    frozen = dict(d=12, k=3, n_values=[50,100,250,500,1000,2000],
                  delta_values=[.05,.10,.20,.35,.50,.75,1.], t_values=[0.,.8],
                  pilot_replications=1000)
    for key, expected in frozen.items():
        if config[key] != expected:
            raise ValueError(f"{key} disagrees with frozen protocol")
    if config['seed_namespaces'] != dict(smoke=0, pilot=1, confirmation_reserved=2):
        raise ValueError("seed namespaces changed")
    return config


def configurations(config):
    return [dict(config_id=f"c{i:03d}", config_index=i, n=n, delta=delta, t=t)
            for i,(n,delta,t) in enumerate(itertools.product(
                config['n_values'], config['delta_values'], config['t_values']))]


def eigenvalues(delta, t):
    if not 0 < delta <= 1 or t not in (0., .8):
        raise ValueError("outside the E3 population family")
    return np.array([6.,4.,2.5+delta/2,2.5-delta/2]+[1+t]*4+[1-t]*4)


def seed_sequence(config, stage, index, purpose):
    # Purpose 0=simulation, 1=joint bootstrap, 2=calibration; no shared draws.
    return np.random.SeedSequence(config['root_seed'], spawn_key=(
        config['seed_namespaces'][stage], index, purpose))


def source_files():
    files = [ROOT/'PCA_stability.py', ROOT/'stability_tests.py']
    for base, pattern in [('src','*.py'), ('tests','*.py'), ('experiments/e3','*.py')]:
        files.extend(sorted((ROOT/base).rglob(pattern)))
    files += [ROOT/'experiments/e3/config.yaml', ROOT/'experiments/e3/E3.md',
              ROOT/'experiments/e3/README.md', ROOT/'requirements-e3.txt']
    return files


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    return {str(p.relative_to(ROOT)):file_hash(p) for p in source_files()}
