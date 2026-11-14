"""Validate -> smoke -> pilot. Confirmation is deliberately unavailable."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'src')]
os.environ.setdefault('MPLCONFIGDIR', '/tmp/e3-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/e3-cache')

import numpy as np
import pandas as pd
import scipy
import matplotlib
import yaml

from experiments.e3.design import (configurations, file_hash, read_config,
                                    seed_sequence, source_hashes)
from experiments.e3.simulation import simulate_configuration


def save_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, default=str)+'\n')


def fingerprint(hashes):
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def validate():
    suite = unittest.TestSuite()
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName('stability_tests'))
    suite.addTests(unittest.defaultTestLoader.discover(str(ROOT/'tests'), pattern='test_*.py'))
    stream = io.StringIO()
    start = time.perf_counter()
    with redirect_stdout(stream):
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    log = stream.getvalue()
    print(log, flush=True)
    if not result.wasSuccessful() or result.skipped:
        raise SystemExit('Tests failed or mandatory tests skipped; simulation blocked')
    print(f'{result.testsRun} tests passed; archived files were not modified.',flush=True)


def audit_records(records, conditions, replications, stage):
    expected = {row['config_id'] for row in conditions}
    if set(records.config_id) != expected or len(records) != len(conditions)*replications:
        raise ValueError('configuration or replication counts do not match')
    if records.duplicated(['config_id','replication']).any():
        raise ValueError('duplicate replication identifiers')
    for row in conditions:
        group = records[records.config_id == row['config_id']]
        if len(group) != replications or set(group.replication) != set(range(replications)):
            raise ValueError('incomplete configuration')
        for key,value in row.items():
            if not (group[key] == value).all():
                raise ValueError(f'inconsistent {key}')
        if not (group.stage == stage).all():
            raise ValueError('incorrect stage')
    if records.failure.any():
        raise ValueError('Numerical failures retained; full analysis blocked pending investigation')


def verify_archive(output):
    """Verify archived numerical evidence without writing or rerunning simulations.

    Original provenance hashes remain historical: removed SVG/report artifacts
    are retained in Git history, not required for the curated archive.
    """
    from experiments.e3.design import eigenvalues
    from PCA_stability import first_order_benchmark
    output=Path(output)
    manifest=json.loads((output/'run_manifest.json').read_text())
    if file_hash(output/'replications.csv')!=manifest['replications_sha256']:
        raise RuntimeError('Original replication file fails SHA-256 verification')
    if file_hash(ROOT/'experiments/e3/E3.md')!=manifest['protocol_sha256']:
        raise RuntimeError('Frozen protocol differs from the original run')
    names=['summary.csv','contrasts.csv','ranking.csv','calibration.csv']
    names += [f'figures/{name}.png' for name in (
        'oracle_signed_error','oracle_absolute_error','oracle_relative_error',
        'relative_median_absolute_error','serious_underprediction',
        'relative_underprediction','prediction_spread')]
    for name in names:
        if file_hash(output/name)!=manifest['analysis_artifact_hashes'][name]:
            raise RuntimeError(f'Numerical artifact is missing or modified: {name}')
    if manifest.get('validation_json_sha256'):
        if file_hash(ROOT/'results/e3/validation.json')!=manifest['validation_json_sha256']:
            raise RuntimeError('Historical validation reference differs')
    if manifest.get('smoke_manifest_sha256'):
        if file_hash(ROOT/'results/e3/smoke/run_manifest.json')!=manifest['smoke_manifest_sha256']:
            raise RuntimeError('Historical smoke reference differs')
    records=pd.read_csv(output/'replications.csv',float_precision='round_trip',keep_default_na=False)
    summary=pd.read_csv(output/'summary.csv',float_precision='round_trip')
    audit_records(records,manifest['conditions'],manifest['replications_per_configuration'],manifest['stage'])
    if len(summary)!=len(manifest['conditions']) or summary.config_id.nunique()!=len(summary):
        raise RuntimeError('Summary grid incomplete or duplicated')
    values=records[['loss','plugin_raw','plugin_capped','sample_gap']].to_numpy()
    if not np.isfinite(values).all() or not records.loss.between(0,1).all():
        raise RuntimeError('Invalid numerical observation')
    if not np.array_equal(records.plugin_capped.to_numpy(),
                          np.minimum(1.,np.maximum(0.,records.plugin_raw.to_numpy()))):
        raise RuntimeError('Capped prediction mismatch')
    def check(actual,expected,label):
        if not np.isclose(actual,expected,rtol=1e-10,atol=1e-12):
            raise RuntimeError('Numerical summary mismatch: '+label)
    for row in summary.itertuples(index=False):
        g=records[records.config_id==row.config_id]
        if len(g)!=row.M: raise RuntimeError('Missing replication for '+row.config_id)
        risk=g.loss.mean()
        benchmark=first_order_benchmark(np.diag(eigenvalues(row.delta,row.t)),
                                        int(row.n),int(manifest['config']['k']))
        check(row.risk,risk,'risk '+row.config_id)
        check(row.oracle,benchmark,'oracle '+row.config_id)
        check(row.median_prediction,g.plugin_raw.median(),'median '+row.config_id)
        check(row.median_absolute_error,np.median(abs(g.plugin_raw-risk)),
              'absolute error '+row.config_id)
        check(row.relative_median_absolute_error,row.median_absolute_error/risk,
              'relative error '+row.config_id)
        check(row.relative_underprediction,np.mean(g.plugin_raw<.5*risk),
              'underprediction '+row.config_id)
        for q in (5,25,50,75,95,99):
            check(getattr(row,f'raw_q{q:02d}'),np.quantile(g.plugin_raw,q/100),
                  f'quantile {q} for {row.config_id}')
    print(f"ARCHIVE VERIFIED (read-only): {len(summary)} configurations, "
          f"{len(records):,} observations, 7 PNG figures and 4 result tables.")
    print('Original source hashes and optional removed artifacts are preserved in Git history.')

def analyze(output):
    from experiments.e3.analysis import summarize, contrasts_and_ranking, calibration
    from experiments.e3.figures import make_figures
    output = Path(output)
    manifest = json.loads((output/'run_manifest.json').read_text())
    if manifest['source_hashes'] != source_hashes():
        raise RuntimeError('Source differs from this run; archived pilot is read-only')
    if manifest['protocol_sha256'] != file_hash(ROOT/'experiments/e3/E3.md'):
        raise RuntimeError('Frozen protocol differs from this run')
    if file_hash(output/'replications.csv') != manifest['replications_sha256']:
        raise RuntimeError('Replication artifact changed')
    records = pd.read_csv(output/'replications.csv', float_precision='round_trip', keep_default_na=False)
    audit_records(records, manifest['conditions'], manifest['replications_per_configuration'], manifest['stage'])
    summary, draws = summarize(records, manifest['config'], manifest['stage'])
    contrasts, ranking = contrasts_and_ranking(summary, draws)
    calibrated = calibration(records, summary, manifest['config'], manifest['stage'])
    for name,frame in [('summary',summary),('contrasts',contrasts),('ranking',ranking),('calibration',calibrated)]:
        frame.to_csv(output/f'{name}.csv', index=False)
    make_figures(summary, output/'figures', manifest['run_id'])
    figures = sorted((output/'figures').glob('*'))
    artifacts = [output/f'{name}.csv' for name in ('summary','contrasts','ranking','calibration')]
    artifacts += figures
    manifest['analysis_artifact_hashes'] = {str(p.relative_to(output)):file_hash(p) for p in artifacts}
    manifest['analysis_complete'] = True
    manifest['verification'] = dict(configurations=len(summary), replications=len(records),
             minimum_count=int(summary.M.min()), maximum_count=int(summary.M.max()),
             numerical_failures=int(records.failure.sum()),
             maximum_decomposition_residual=float(summary.decomposition_residual.abs().max()))
    save_json(output/'run_manifest.json', manifest)
    print(f'Analysis complete: {len(summary)} configurations, {len(records):,} matched records -> {output}', flush=True)


def simulate(stage, output):
    """Simulate only into a fresh directory; preserve the original pilot."""
    if stage not in ('smoke','pilot'):
        raise ValueError('only exploratory stages are implemented')
    validate()
    config=read_config()
    all_conditions=configurations(config)
    conditions=all_conditions if stage=='pilot' else [
        c for c in all_conditions if (c['n'],c['delta']) in ((50,.05),(2000,1.))]
    replications=config[stage+'_replications']
    if stage=='pilot':
        for condition in [c for c in all_conditions
                          if (c['n'],c['delta']) in ((50,.05),(2000,1.))]:
            first=simulate_configuration(config,condition,'smoke',config['smoke_replications'])
            repeat=simulate_configuration(config,condition,'smoke',config['smoke_replications'])
            if not first.equals(repeat) or first.failure.any():
                raise RuntimeError('Smoke reproducibility check failed; pilot blocked')
        print('Four smoke configurations passed.',flush=True)
    output = Path(output)
    if output.exists():
        raise FileExistsError(f'Existing output: {output}; choose a fresh directory')
    output.mkdir(parents=True,exist_ok=False)
    start = time.perf_counter()
    frames = []
    reproducible = True
    for condition in conditions:
        frame = simulate_configuration(config, condition, stage, replications)
        frames.append(frame)
        if stage == 'smoke':
            repeat = simulate_configuration(config, condition, stage, replications)
            reproducible &= frame.equals(repeat)
        print(f"{stage} {condition['config_id']}: n={condition['n']}, delta={condition['delta']}, "
              f"t={condition['t']}, M={len(frame)}, failures={frame.failure.sum()}", flush=True)
    records = pd.concat(frames, ignore_index=True)
    records.to_csv(output/'replications.csv', index=False)
    # Round-trip through the actual persisted artifact, preserving float precision.
    reread = pd.read_csv(output/'replications.csv', float_precision='round_trip', keep_default_na=False)
    pd.testing.assert_frame_equal(records, reread, check_dtype=False, check_exact=True)
    hashes = source_hashes()
    revision = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    status = subprocess.check_output(['git','status','--porcelain'], cwd=ROOT, text=True)
    run_id = f"{stage}-{fingerprint(hashes)[:12]}-{config['root_seed']}"
    manifest = dict(stage=stage, run_id=run_id, config=config, conditions=conditions,
        replications_per_configuration=replications, total_replications=len(records),
        source_hashes=hashes, source_fingerprint=fingerprint(hashes),
        git_base_commit=revision, git_status_porcelain=status,
        protocol_sha256=file_hash(ROOT/'experiments/e3/E3.md'),
        tests_passed=True, smoke_checked=(stage=='pilot'),
        replications_sha256=file_hash(output/'replications.csv'),
        versions=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                      pandas=pd.__version__, matplotlib=matplotlib.__version__, pyyaml=yaml.__version__),
        environment=dict(platform=platform.platform(), machine=platform.machine(),
                         cpu_count=os.cpu_count(), executable=sys.executable,
                         thread_settings={k:os.environ.get(k) for k in
                             ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
                         numpy_blas=str(np.__config__.CONFIG)),
        utc=datetime.now(timezone.utc).isoformat(), simulation_seconds=time.perf_counter()-start,
        numerical_failures=int(records.failure.sum()),
        seed_strategy='PCG64(SeedSequence(root_seed, spawn_key=(stage,config_index,purpose))); independent configurations',
        streams={c['config_id']:{p:list(seed_sequence(config,stage,c['config_index'],j).spawn_key)
                                for p,j in [('simulation',0),('joint_bootstrap',1),
                                            ('calibration_if_first_in_stratum',2)]} for c in conditions},
        analysis_complete=False)
    save_json(output/'run_manifest.json', manifest)
    audit_records(reread, conditions, replications, stage)
    if stage == 'smoke' and not reproducible:
        raise RuntimeError('Smoke seed reproducibility failed')
    analyze(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate','smoke','pilot','analyze','verify'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.command == 'validate':
        if args.output:
            parser.error('validate does not use --output')
        validate()
    elif args.command == 'verify':
        verify_archive(args.output or ROOT/'results/e3/pilot')
    elif args.command == 'analyze':
        if not args.output:
            parser.error('analyze requires --output (do not edit the archived pilot)')
        analyze(args.output)
    else:
        if not args.output:
            parser.error('simulations require a new --output directory')
        simulate(args.command,args.output)
