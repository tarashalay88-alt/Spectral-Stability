"""Independent postflight audit; run from the Spectral-Stability repository root."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

root=Path.cwd()
sys.path[:0]=[str(root),str(root/'src')]
from experiments.e3.design import source_hashes, eigenvalues, configurations, read_config, file_hash
from experiments.e3.simulation import simulate_configuration

base=root/'results/e3'
p=base/'pilot'
m=json.loads((p/'run_manifest.json').read_text())
s=pd.read_csv(p/'summary.csv',float_precision='round_trip')
r=pd.read_csv(p/'replications.csv',float_precision='round_trip',keep_default_na=False)
checks={}
checks['source_hashes_current']=m['source_hashes']==source_hashes()
checks['protocol_unchanged']=subprocess.check_output(['git','show',m['git_base_commit']+':experiments/e3/E3.md'])==(root/'experiments/e3/E3.md').read_bytes()
checks['all_artifact_hashes']=all(file_hash(p/name)==digest for name,digest in m['analysis_artifact_hashes'].items())
checks['replications_hash']=file_hash(p/'replications.csv')==m['replications_sha256']
checks['validation_reference']=file_hash(base/'validation.json')==m['validation_json_sha256']
checks['smoke_reference']=file_hash(base/'smoke/run_manifest.json')==m['smoke_manifest_sha256']
checks['counts_and_grid']=len(s)==84 and len(r)==84000 and (r.groupby('config_id').size()==1000).all() and len(configurations(read_config()))==84
checks['unique_ids']=not r.duplicated(['config_id','replication']).any()
checks['finite_all_records']=np.isfinite(r[['loss','plugin_raw','plugin_capped','sample_gap']]).all().all()
checks['zero_failures']=not r.failure.any()
checks['actual_losses_in_unit_interval']=r.loss.between(0,1).all()
checks['caps_exact']=np.array_equal(r.plugin_capped,np.minimum(1,np.maximum(0,r.plugin_raw)))
checks['extremes_retained']=int((r.plugin_raw>1).sum())==3032 and r.plugin_raw.max()>8000
checks['unique_simulation_keys']=r.groupby('config_id').seed_spawn_key.first().nunique()==84
for row in s.itertuples():
    g=r[r.config_id==row.config_id]
    assert len(g)==1000 and set(g.replication)==set(range(1000))
    assert np.all(g.n==row.n) and np.all(g.delta==row.delta) and np.all(g.t==row.t)
    risk=float(math.fsum(g.loss)/1000)
    se=float(np.sqrt(np.sum((g.loss-risk)**2)/999)/np.sqrt(1000))
    assert np.isclose(row.risk,risk,rtol=1e-14)
    assert np.isclose(row.risk_mcse,se,rtol=1e-14)
    lam=eigenvalues(row.delta,row.t)
    oracle=sum(lam[i]*lam[j]/(lam[i]-lam[j])**2 for i in range(3) for j in range(3,12))/(3*(row.n-1))
    assert np.isclose(row.oracle,oracle,rtol=1e-13)
    # Independent closed-form accepted-risk interval, rather than piecewise evaluation.
    accept_lo=min(oracle-.002,oracle/1.2)
    accept_hi=max(oracle+.002,oracle/.8)
    label=('within' if row.risk_low>=accept_lo and row.risk_high<=accept_hi else
           'outside' if row.risk_high<accept_lo or row.risk_low>accept_hi else 'uncertain')
    assert label==row.oracle_label
    assert np.isclose(row.median_absolute_error,np.median(abs(g.plugin_raw-risk)),rtol=1e-13)
    assert np.isclose(row.relative_median_absolute_error,row.median_absolute_error/risk,rtol=1e-13)
    assert np.isclose(row.capped_mae,np.abs(g.plugin_capped-risk).mean(),rtol=1e-13)
    assert np.isclose(row.capped_rmse,np.sqrt(np.mean((g.plugin_capped-risk)**2)),rtol=1e-13)
    for q in (5,25,50,75,95,99):
        assert np.isclose(getattr(row,f'raw_q{q:02d}'),np.quantile(g.plugin_raw,q/100),rtol=1e-13)
    assert np.isclose(row.serious_underprediction,np.mean(g.plugin_raw<=.05))
    assert np.isclose(row.relative_underprediction,np.mean(g.plugin_raw<.5*risk))
    assert abs(row.decomposition_residual)<1e-12
checks['independent_summary_formulas_and_oracle_labels']=True
# Recover complete original pilot draws for two configurations, including the largest extreme.
config=read_config()
for index in (0,43):
    generated=simulate_configuration(config,configurations(config)[index],'pilot',1000)
    saved=r[r.config_index==index].reset_index(drop=True)
    pd.testing.assert_frame_equal(generated,saved,check_dtype=False,check_exact=True)
checks['selected_pilot_streams_exactly_recovered']=True
c=pd.read_csv(p/'calibration.csv')
checks['secondary_complete']=(len(pd.read_csv(p/'contrasts.csv'))==42 and len(pd.read_csv(p/'ranking.csv'))==12 and len(c)==60 and (c.groupby(['n','t'])['count'].sum()==7000).all())
checks['figures_complete']=len(list((p/'figures').glob('*.png')))==7 and len(list((p/'figures').glob('*.svg')))==7
previous=json.loads((p/'regeneration_reference.json').read_text())
checks['regenerated_tables_identical']=all(previous[name]==m['analysis_artifact_hashes'][name] for name in ['summary.csv','contrasts.csv','ranking.csv','calibration.csv','pilot_report.md'])
checks={k:bool(v) for k,v in checks.items()}
assert all(checks.values()),checks
output=dict(passed=True,checks=checks,run_id=m['run_id'],replications=84000,configurations=84,replications_each=1000,raw_exceeds_one=3032,maximum_raw=float(r.plugin_raw.max()),maximum_decomposition_residual=float(s.decomposition_residual.abs().max()))
(p/'audit.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
