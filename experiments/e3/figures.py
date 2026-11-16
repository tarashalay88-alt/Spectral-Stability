"""Plots only simulated cells; metadata links every figure to its manifest."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/e3-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/e3-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def serious_underprediction_color_limits(summary):
    """Scale only on eligible cells, never on masked/ineligible probabilities."""
    eligible = summary.loc[
        summary['serious_eligibility'] == 'established', 'serious_underprediction'
    ].to_numpy(dtype=float)
    eligible = eligible[np.isfinite(eligible)]
    if not len(eligible):
        return 0.0, 1.0
    # Round up to a readable tenth; exclude values that will not be plotted.
    return 0.0, max(0.1, np.ceil(float(eligible.max()) * 10) / 10)


def make_figures(summary, output, run_id):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    ns, deltas = sorted(summary.n.unique()), sorted(summary.delta.unique())
    footer = f'Exploratory E3 {run_id}; d=12, k=3. Seeds/code/software: ../run_manifest.json'
    maps = [
        ('oracle_signed_error', 'Oracle signed error B − R̂', False),
        ('oracle_absolute_error', 'Oracle absolute error |B − R̂|', True),
        ('oracle_relative_error', 'Oracle relative error |B − R̂| / R̂', True),
        ('relative_median_absolute_error', 'Median raw plug-in absolute error / R̂', True),
        ('serious_underprediction', r'Serious underprediction $P(\\widehat{B}\\leq 0.05)$; risk CI above 0.10', False),
        ('relative_underprediction', 'Relative underprediction P(B̂ < 0.5 R̂)', False),
    ]
    for metric, title, logarithmic in maps:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
        observed = summary[metric].to_numpy()
        if logarithmic:
            positives = observed[observed > 0]
            norm = LogNorm(vmin=max(1e-6, positives.min()), vmax=max(positives.max(), positives.min()*1.01))
        else:
            norm = None
        kwargs = dict(norm=norm) if logarithmic else dict(vmin=min(0., observed.min()), vmax=max(observed.max(), .001))
        if metric == 'oracle_signed_error':
            bound = max(abs(observed.min()), abs(observed.max()))
            kwargs = dict(vmin=-bound, vmax=bound)
        if metric == 'serious_underprediction':
            # Masked cells can have probabilities close to 1.0; they MUST NOT
            # determine the displayed color range for eligible cells.
            lower, upper = serious_underprediction_color_limits(summary)
            kwargs = dict(vmin=lower, vmax=upper)
        for ax, t in zip(axes, (0., .8)):
            subset = summary[summary.t == t]
            cells = np.full((len(deltas), len(ns)), np.nan)
            for row in subset.to_dict('records'):
                i,j = deltas.index(row['delta']), ns.index(row['n'])
                if metric != 'serious_underprediction' or row['serious_eligibility'] == 'established':
                    cells[i,j] = row[metric]
                else:
                    ax.text(j,i,'?' if row['serious_eligibility']=='uncertain' else '—',
                            ha='center',va='center',fontsize=8)
            im = ax.imshow(np.ma.masked_invalid(cells), origin='lower', aspect='auto',
                           cmap=('coolwarm' if metric == 'oracle_signed_error' else
                                 plt.get_cmap('viridis').copy().with_extremes(bad='#f1f1f1')),
                           **kwargs)
            ax.set(xticks=range(len(ns)), xticklabels=ns, yticks=range(len(deltas)),
                   yticklabels=[f'{d:.2f}' for d in deltas], xlabel='Sample size n',
                   ylabel='Population boundary eigengap δ', title=f'Background t={t:g}')
        cb_label = ('Probability (eligible configurations only)'
                    if metric == 'serious_underprediction'
                    else 'Log color scale' if logarithmic else 'Value')
        fig.colorbar(im, ax=axes, shrink=.85, label=cb_label)
        note = ('\nEligibility: lower 95% CI for R > 0.10; '
                '— = below threshold, ? = uncertain'
                if metric == 'serious_underprediction' else '')
        fig.suptitle(title+'\n'+footer+note, fontsize=10)
        fig.savefig(output/f'{metric}.png',dpi=160)
        plt.close(fig)
    fig, axes = plt.subplots(len(ns), 2, figsize=(11, 2.6*len(ns)), squeeze=False,
                             layout='constrained')
    for i,n in enumerate(ns):
        for j,t in enumerate((0.,.8)):
            ax = axes[i,j]
            subset = summary[(summary.n==n)&(summary.t==t)].sort_values('delta')
            x = subset.delta.to_numpy()
            ax.fill_between(x, subset.raw_q05, subset.raw_q95, alpha=.18, label='Raw 5–95%')
            ax.fill_between(x, subset.raw_q25, subset.raw_q75, alpha=.35, label='Raw 25–75%')
            ax.vlines(x, subset.raw_q05, subset.raw_q95, alpha=.35, linewidth=1)
            ax.vlines(x, subset.raw_q25, subset.raw_q75, alpha=.5, linewidth=3)
            ax.plot(x, subset.raw_q50, 'o-', label='Raw median')
            ax.plot(x, subset.raw_q99, 'x:', label='Raw 99%')
            ax.plot(x, subset.risk, 's--', color='black', label='Risk R̂')
            ax.set(yscale='log', xlabel='Population boundary eigengap δ',
                   ylabel='Raw B̂ or risk R̂', title=f'n={n}, t={t:g}')
            ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7, ncols=2)
    fig.suptitle('Raw plug-in prediction spread (log y scale)\n'+footer, fontsize=10)
    fig.savefig(output/'prediction_spread.png',dpi=160)
    plt.close(fig)
