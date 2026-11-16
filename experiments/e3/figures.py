"""E3 figures from saved summaries; no simulations or estimates are changed.

Run in PyCharm and choose summary.csv, or use:
    python figures.py --run-dir PATH --no-show
PNG and SVG figures are saved beside this script unless --output is supplied.
Requires numpy, pandas, matplotlib.
"""
from pathlib import Path
import hashlib
import json
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm, SymLogNorm


def serious_underprediction_color_limits(summary):
    values = summary.loc[summary.serious_eligibility == 'established',
                         'serious_underprediction'].to_numpy(float)
    values = values[np.isfinite(values)]
    return 0., max(.1, np.ceil(values.max()*10)/10) if len(values) else 1.


def oracle_sign_classes(summary):
    """Pointwise 95% interval decision: -1 negative, 0 unresolved, +1 positive."""
    return np.where(summary.oracle > summary.risk_high, 1,
                    np.where(summary.oracle < summary.risk_low, -1, 0))


def _validate(summary):
    keys = ['n', 'delta', 't']
    if summary.duplicated(keys).any():
        raise ValueError('Duplicate configurations in summary.csv')
    if summary.empty or not set(summary.t.unique()).issubset({0., .8}):
        raise ValueError('This E3 layout requires background t=0 and/or t=0.8')
    required = keys + ['risk', 'risk_low', 'risk_high', 'oracle',
        'oracle_signed_error', 'oracle_absolute_error', 'oracle_relative_error',
        'relative_median_absolute_error', 'relative_underprediction',
        'serious_underprediction', 'serious_eligibility']
    for q in (5,25,50,75,95,99):
        required += [f'raw_q{q:02d}', f'raw_q{q:02d}_low', f'raw_q{q:02d}_high']
    for name in ('relative_median_absolute_error','relative_underprediction',
                 'serious_underprediction'):
        required += [name+'_low', name+'_high']
    missing = set(required)-set(summary.columns)
    if missing:
        raise ValueError('Missing columns: '+', '.join(sorted(missing)))
    if not np.isfinite(summary[keys+['risk','oracle']].to_numpy(float)).all():
        raise ValueError('Nonfinite configuration, risk or oracle value')
    if (summary.risk <= 0).any():
        raise ValueError('Relative metrics require positive risk')
    if 'failure_count' in summary and (summary.failure_count != 0).any():
        raise ValueError('Resolve recorded numerical failures before plotting')
    if ((summary.risk_low > summary.risk) | (summary.risk_high < summary.risk)).any():
        raise ValueError('Invalid risk intervals')
    q = summary[[f'raw_q{x:02d}' for x in (5,25,50,75,95,99)]].to_numpy(float)
    if not np.isfinite(q).all() or (q <= 0).any() or (np.diff(q,axis=1)<0).any():
        raise ValueError('Log plots require finite positive ordered quantiles')
    expected = np.where(summary.risk_low > .10, 'established',
                        np.where(summary.risk_high < .10, 'below_threshold', 'uncertain'))
    if not np.array_equal(expected, summary.serious_eligibility):
        raise ValueError('Serious-underprediction eligibility differs from risk intervals')


def _save(fig, output, name, footer, note, show):
    fig.text(.5,.012, note+'\n'+footer, ha='center', va='bottom', fontsize=8)
    for extension in ('png','svg'):
        fig.savefig(output/f'{name}.{extension}', dpi=180, bbox_inches='tight')
    if not show:
        plt.close(fig)


def _heatmap(summary, output, name, title, footer, show, ns, deltas, ts):
    sign = name == 'oracle_signed_error_sign'
    source = 'oracle_signed_error' if sign else name
    observed = summary[source].to_numpy(float)
    finite = observed[np.isfinite(observed)]
    positive = finite[finite>0]
    cmap = plt.get_cmap('viridis').with_extremes(bad='#eeeeee')
    norm = None
    cb_label = title.split('$')[0].strip()
    note = 'Discrete tested grid; equal cell spacing is categorical.'
    if sign:
        cmap = ListedColormap(['#2166ac','#eeeeee','#b2182b'])
        norm = BoundaryNorm([-1.5,-.5,.5,1.5],3)
        note += '\nGrey: sign unresolved by pointwise 95% Monte Carlo risk CI; not evidence of zero error.'
    elif name == 'oracle_signed_error':
        bound = max(float(np.max(abs(finite))), 1e-12)
        # Resolve small errors while retaining all extremes; report the scale.
        threshold = max(float(np.median(summary.risk_mcse)), 1e-8) if 'risk_mcse' in summary else 1e-4
        norm = SymLogNorm(linthresh=threshold, vmin=-bound, vmax=bound, base=10)
        cmap = plt.get_cmap('coolwarm').with_extremes(bad='#eeeeee')
        cb_label = f'Signed error (symmetric log; linear within ±{threshold:.2g})'
        note += '\nColours show point estimates; use the sign map for uncertainty.'
    elif name in ('oracle_absolute_error','oracle_relative_error','relative_median_absolute_error'):
        if not len(positive):
            norm = plt.Normalize(0,1)
        else:
            norm = LogNorm(positive.min(), max(positive.max(),positive.min()*1.01))
        cb_label = 'Error (logarithmic colour scale)'
        note += '\nHatcing: pointwise 95% interval includes B (oracle maps) or the 20% relative-error reference (plug-in map).'
        note = note.replace('Hatcing','Hatching')
    else:
        upper = serious_underprediction_color_limits(summary)[1] if name == 'serious_underprediction' else 1.
        norm = plt.Normalize(0,upper)
        cb_label = 'Estimated probability' + (' (eligible configurations)' if name == 'serious_underprediction' else '')
        if name == 'serious_underprediction':
            note += '\nEligibility: lower 95% risk CI > 0.10; —: upper CI < 0.10; ?: unresolved. Coloured cells are eligible only.'
    fig, axes = plt.subplots(1,len(ts),figsize=(12,5.6),squeeze=False)
    fig.subplots_adjust(left=.075,right=.88,bottom=.25,top=.83,wspace=.25)
    for ax,t in zip(axes[0],ts):
        cells = np.full((len(deltas),len(ns)),np.nan)
        group = summary[summary.t == t].copy()
        group['_sign'] = oracle_sign_classes(group)
        for row in group.to_dict('records'):
            i,j = deltas.index(row['delta']),ns.index(row['n'])
            eligible = name != 'serious_underprediction' or row['serious_eligibility']=='established'
            if not eligible:
                ax.text(j,i,'?' if row['serious_eligibility']=='uncertain' else '—',ha='center',va='center')
                continue
            cells[i,j] = row['_sign'] if sign else row[source]
            hatch = (name in ('oracle_absolute_error','oracle_relative_error') and row['risk_low'] <= row['oracle'] <= row['risk_high'])
            hatch |= (name == 'relative_median_absolute_error' and row[source+'_low'] <= .20 <= row[source+'_high'])
            if hatch:
                from matplotlib.patches import Rectangle
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,fill=False,hatch='///',edgecolor='#666666',linewidth=0,zorder=3))
            if not sign and (name.endswith('underprediction') or name=='relative_median_absolute_error'):
                value = row[source]
                text = f'{value:.3f}\n[{row[source+"_low"]:.3f}, {row[source+"_high"]:.3f}]'
                rgba = cmap(norm(value)); luminance = .2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                ax.text(j,i,text,ha='center',va='center',fontsize=6.6,color='black' if luminance>.52 else 'white',zorder=4)
        im = ax.imshow(np.ma.masked_invalid(cells),origin='lower',aspect='auto',cmap=cmap,norm=norm)
        ax.set(xticks=range(len(ns)),xticklabels=ns,yticks=range(len(deltas)),yticklabels=[f'{d:.2f}' for d in deltas],xlabel='Sample size n',ylabel='Population boundary eigengap δ',title=f'Background t={t:g}')
    cb = fig.colorbar(im,cax=fig.add_axes([.905,.27,.016,.52]),label=cb_label)
    if sign:
        cb.set_ticks([-1,0,1],labels=['Negative','Unresolved','Positive'])
        cb.set_label('Sign supported by pointwise 95% risk CI')
    if name.endswith('underprediction') or name == 'relative_median_absolute_error':
        note += '\nCell annotations: estimate [pointwise 95% CI]; intervals are not simultaneous.'
    fig.suptitle(title,fontsize=15,y=.96)
    _save(fig,output,name,footer,note,show)


def _plot_spread(summary, ns, footer, path, show=False):
    ts = sorted(summary.t.unique())
    fig,axes = plt.subplots(len(ns),len(ts),figsize=(12,3*len(ns)+1.3),squeeze=False)
    fig.subplots_adjust(left=.08,right=.98,bottom=.12 if len(ns)>3 else .18,top=.88 if len(ns)<=3 else .94,hspace=.55,wspace=.25)
    for i,n in enumerate(ns):
        for j,t in enumerate(ts):
            ax=axes[i,j]
            sub=summary[(summary.n==n)&(summary.t==t)].sort_values('delta')
            x=sub.delta.to_numpy()
            ax.fill_between(x,sub.raw_q05,sub.raw_q95,color='#4393c3',alpha=.22,label='Sample distribution: 5–95%')
            ax.fill_between(x,sub.raw_q25,sub.raw_q75,color='#f4a259',alpha=.5,label='Sample distribution: 25–75%')
            ax.plot(x,sub.raw_q50,'o-',color='#2166ac',label='Raw median')
            ax.errorbar(x,sub.raw_q99,yerr=np.array([sub.raw_q99-sub.raw_q99_low,sub.raw_q99_high-sub.raw_q99]),fmt='x:',color='#e87817',capsize=2,label='Raw 99th percentile + 95% CI')
            ax.errorbar(x,sub.risk,yerr=np.array([sub.risk-sub.risk_low,sub.risk_high-sub.risk]),fmt='s--',color='black',capsize=3,label=r'Risk $\widehat R$ + 95% MC CI')
            ax.set(yscale='log',xlabel='Population boundary eigengap δ',ylabel=r'Raw $\widehat B$ or risk $\widehat R$',title=f'n={n}, t={t:g}')
            ax.grid(alpha=.2)
        # Shared scale within each n: direct comparison between backgrounds.
        bounds=[ax.get_ylim() for ax in axes[i]]
        for ax in axes[i]: ax.set_ylim(min(v[0] for v in bounds),max(v[1] for v in bounds))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.968 if len(ns)>3 else .955),ncol=3,fontsize=8,frameon=False)
    fig.suptitle('Raw plug-in prediction spread (log scale)',y=.995,fontsize=16)
    _save(fig,Path(path).parent,Path(path).stem,footer,
          'Shading: variation across samples, not confidence bands. Error bars: pointwise Monte Carlo uncertainty.\nY limits shared across backgrounds within each row; may differ between sample sizes.',show)


def make_figures(summary, output, run_id=None, show=False, metadata=None):
    _validate(summary)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    metadata=metadata or {}
    candidate=output.parent/'run_manifest.json'
    if not metadata and candidate.is_file():
        metadata=json.loads(candidate.read_text(encoding='utf-8'))
        metadata['_manifest_path']=str(candidate.resolve())
    config=metadata.get('config',{})
    identity=metadata.get('run_id',run_id) or 'unverified run identity'
    dims=', '.join(f'{key}={config[key]}' for key in ('d','k') if key in config)
    m = ', '.join(map(str,sorted(summary.M.unique()))) if 'M' in summary else 'unknown'
    footer=f'Exploratory E3 · {identity}'+(f' · {dims}' if dims else '')+f' · M={m} per cell'
    manifest_path=metadata.get('_manifest_path')
    reference=os.path.relpath(manifest_path,output.resolve()) if manifest_path else 'run manifest not supplied'
    footer+='\nSimulation parameters, seeds and software: '+reference
    maps=[('oracle_signed_error',r'Oracle signed error $B-\widehat R$'),
          ('oracle_signed_error_sign',r'Oracle error sign with Monte Carlo uncertainty'),
          ('oracle_absolute_error',r'Oracle absolute error $|B-\widehat R|$'),
          ('oracle_relative_error',r'Oracle relative error $|B-\widehat R|/\widehat R$'),
          ('relative_median_absolute_error',r'Median raw plug-in absolute error / $\widehat R$'),
          ('serious_underprediction',r'Serious underprediction $P(\widehat B\leq0.05)$, eligible risks only'),
          ('relative_underprediction',r'Relative underprediction $P(\widehat B<0.5 R)$ (estimated using $\widehat R$)')]
    ns=sorted(summary.n.unique());deltas=sorted(summary.delta.unique());ts=sorted(summary.t.unique())
    for name,title in maps: _heatmap(summary,output,name,title,footer,show,ns,deltas,ts)
    main=[n for n in (100,500,2000) if n in ns]
    _plot_spread(summary,main or ns,footer,output/'prediction_spread.png',show)
    _plot_spread(summary,ns,footer,output/'prediction_spread_supplement.png',show)



if __name__ == '__main__':
    import argparse
    import pandas as pd
    parser=argparse.ArgumentParser(description=__doc__)
    script=Path(__file__).resolve()
    default_run=next((p/'results/e3/pilot' for p in script.parents if (p/'results/e3/pilot/summary.csv').is_file()),None)
    # Bundled results are convenient when running the downloaded ZIP.
    if (script.parent/'summary.csv').is_file(): default_run=script.parent
    parser.add_argument('--run-dir',type=Path,default=default_run)
    parser.add_argument('--output',type=Path,default=script.parent)
    parser.add_argument('--no-show',action='store_true')
    args=parser.parse_args()
    if args.run_dir is None:
        from tkinter import Tk,filedialog
        root=Tk();root.withdraw()
        try: selected=filedialog.askopenfilename(title='Select saved E3 summary.csv',filetypes=[('CSV files','*.csv')])
        finally: root.destroy()
        if not selected: parser.exit(message='No summary file selected.\n')
        summary_path=Path(selected)
    else: summary_path=args.run_dir.resolve()/'summary.csv'
    summary=pd.read_csv(summary_path,float_precision='round_trip')
    manifest=summary_path.parent/'run_manifest.json'
    metadata=json.loads(manifest.read_text(encoding='utf-8')) if manifest.is_file() else {}
    digest=hashlib.sha256(summary_path.read_bytes()).hexdigest()
    expected=metadata.get('analysis_artifact_hashes',{}).get('summary.csv')
    if expected and digest!=expected:
        raise ValueError('summary.csv hash differs from its run manifest; choose the matching files')
    if not manifest.is_file(): warnings.warn('No run_manifest.json found; run identity and dimensions are unavailable')
    metadata['_manifest_path']=str(manifest.resolve()) if manifest.is_file() else None
    metadata['_source']={'summary_file':summary_path.name,'summary_sha256':digest,
                         'manifest_file':manifest.name if manifest.is_file() else None,
                         'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest() if manifest.is_file() else None,
                         'summary_hash_matches_manifest':digest==expected if expected else None}
    make_figures(summary,args.output,show=not args.no_show,metadata=metadata)
    print(f'Figures saved beside figures.py: {args.output.resolve()}')
    if not args.no_show: plt.show()
