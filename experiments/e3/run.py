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
    output = ROOT/'results/e3'
    output.mkdir(parents=True, exist_ok=True)
    suite = unittest.TestSuite()
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName('stability_tests'))
    suite.addTests(unittest.defaultTestLoader.discover(str(ROOT/'tests'), pattern='test_*.py'))
    stream = io.StringIO()
    start = time.perf_counter()
    with redirect_stdout(stream):
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    log = stream.getvalue()
    print(log, flush=True)
    (output/'validation.log').write_text(log)
    hashes = source_hashes()
    save_json(output/'validation.json', dict(passed=result.wasSuccessful(),
              tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
              skipped=len(result.skipped), source_hashes=hashes, source_fingerprint=fingerprint(hashes),
              utc=datetime.now(timezone.utc).isoformat(), seconds=time.perf_counter()-start,
              log_sha256=file_hash(output/'validation.log')))
    if not result.wasSuccessful() or result.skipped:
        raise SystemExit('Validation failed or mandatory tests skipped; pilot blocked')


def require_validation():
    path = ROOT/'results/e3/validation.json'
    if not path.exists():
        raise RuntimeError('Run validate first')
    gate = json.loads(path.read_text())
    if not gate['passed'] or gate['skipped'] or gate['source_hashes'] != source_hashes():
        raise RuntimeError('Validation missing, failed, skipped, or stale; run validate')
    if file_hash(path.with_suffix('.log')) != gate['log_sha256']:
        raise RuntimeError('Validation log changed')
    return path


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


def analyze(output):
    from experiments.e3.analysis import summarize, contrasts_and_ranking, calibration
    from experiments.e3.figures import make_figures
    from experiments.e3.report import make_report
    output = Path(output)
    manifest = json.loads((output/'run_manifest.json').read_text())
    if manifest['source_hashes'] != source_hashes():
        raise RuntimeError('Source differs from simulation; review changes before regenerating')
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
    make_report(summary, records, manifest, output/'pilot_report.md')
    figures = sorted((output/'figures').glob('*'))
    save_json(output/'figure_index.json', dict(run_id=manifest['run_id'],
              source_fingerprint=manifest['source_fingerprint'],
              manifest='run_manifest.json', figures={p.name:file_hash(p) for p in figures}))
    artifacts = [output/f'{name}.csv' for name in ('summary','contrasts','ranking','calibration')]
    artifacts += figures+[output/'pilot_report.md', output/'figure_index.json']
    manifest['analysis_artifact_hashes'] = {str(p.relative_to(output)):file_hash(p) for p in artifacts}
    manifest['analysis_complete'] = True
    manifest['verification'] = dict(configurations=len(summary), replications=len(records),
             minimum_count=int(summary.M.min()), maximum_count=int(summary.M.max()),
             numerical_failures=int(records.failure.sum()),
             maximum_decomposition_residual=float(summary.decomposition_residual.abs().max()))
    save_json(output/'run_manifest.json', manifest)
    print(f'Analysis complete: {len(summary)} configurations, {len(records):,} matched records -> {output}', flush=True)


def simulate(stage, output):
    validation_path = require_validation()
    config = read_config()
    all_conditions = configurations(config)
    conditions = all_conditions if stage == 'pilot' else [c for c in all_conditions
        if (c['n'],c['delta']) in ((50,.05),(2000,1.))]
    replications = config[stage+'_replications']
    smoke_path = ROOT/'results/e3/smoke/run_manifest.json'
    if stage == 'pilot':
        if not smoke_path.exists():
            raise RuntimeError('Run default smoke before pilot')
        smoke = json.loads(smoke_path.read_text())
        if (not smoke.get('smoke_passed') or smoke['source_hashes'] != source_hashes()
                or not smoke.get('analysis_complete')):
            raise RuntimeError('Smoke failed, incomplete, or stale')
        for name,digest in smoke['analysis_artifact_hashes'].items():
            if file_hash(smoke_path.parent/name) != digest:
                raise RuntimeError('Smoke analysis artifact changed')
        if file_hash(smoke_path.parent/'replications.csv') != smoke['replications_sha256']:
            raise RuntimeError('Smoke replication artifact changed')
    output = Path(output)
    if (output/'run_manifest.json').exists():
        raise FileExistsError(f'Existing run at {output}; use a fresh --output directory')
    output.mkdir(parents=True, exist_ok=True)
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
        validation_json_sha256=file_hash(validation_path), validation_path='results/e3/validation.json',
        smoke_manifest_sha256=file_hash(smoke_path) if stage=='pilot' else None,
        smoke_path='results/e3/smoke/run_manifest.json' if stage=='pilot' else None,
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
        analysis_complete=False, smoke_passed=False)
    save_json(output/'run_manifest.json', manifest)
    audit_records(reread, conditions, replications, stage)
    if stage == 'smoke' and not reproducible:
        raise RuntimeError('Smoke seed reproducibility failed')
    analyze(output)
    if stage == 'smoke':
        manifest = json.loads((output/'run_manifest.json').read_text())
        manifest['smoke_passed'] = True
        manifest['smoke_checks'] = dict(repeated_seeds_identical=reproducible,
              csv_round_trip_exact=True, full_pipeline_passed=True,
              raw_predictions_above_one_retained=int((records.plugin_raw>1).sum()))
        save_json(output/'run_manifest.json', manifest)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate','smoke','pilot','analyze'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.command == 'validate':
        if args.output:
            parser.error('validation uses results/e3/validation.json')
        validate()
    elif args.command == 'analyze':
        analyze(args.output or ROOT/'results/e3/pilot')
    else:
        simulate(args.command, args.output or ROOT/'results/e3'/args.command)
