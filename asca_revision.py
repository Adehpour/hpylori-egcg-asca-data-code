"""Independent ASCA reimplementation for the Darvish Nejad revision.

Run: python asca_revision.py --input asca_input.csv --out results --permutations 10000
Input has one row per strain, treatment and replicate and the 12 named endpoints.
Published permutation P values are not used in the computation.
"""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STRAINS = ['26695', 'CLR-resistant']
TREATMENTS = ['Untreated', 'Blank NP', 'Free EGCG', 'EGCG-NP', 'Clarithromycin']
GENES = ['cagA', 'vacA', 'ureA', 'flaA', 'babA', 'sabA']
FIELDS = ['Urease_pct', 'Biofilm_pct', *[g + '_dCt' for g in GENES],
          'Associated_log10_CFU_well', 'Phospho_total_CagA', 'Nuc_cyto_p65', 'IL8_pg_mL']
LABELS = ['Urease', 'Biofilm', *GENES, 'Cell-assoc. CFU', 'pCagA/CagA', 'Nuc/cyto p65', 'log₁₀ IL-8']
COLORS = ['#465366', '#8b6c4a', '#32799b', '#137966', '#ae574d']

def data_cube(path):
    source = pd.read_csv(path, dtype={'strain': str})
    keys = ['strain', 'treatment', 'replicate']
    expected = pd.MultiIndex.from_product([STRAINS, TREATMENTS, range(1, 7)], names=keys)
    source = source.set_index(keys)
    if not source.index.is_unique or len(source) != 60 or set(source.index) != set(expected):
        raise ValueError('Expected 60 unique rows covering the specified balanced design.')
    data = source.reindex(expected)[FIELDS].to_numpy(dtype=float)
    if not np.isfinite(data).all() or (data[:, -1] <= 0).any():
        raise ValueError('An endpoint is missing, non-finite or invalid for log transformation.')
    data[:, 2:8] *= -1
    data[:, -1] = np.log10(data[:, -1])
    sd = data.std(axis=0, ddof=1)
    if (sd <= 0).any():
        raise ValueError('All endpoints must vary.')
    return ((data - data.mean(axis=0)) / sd).reshape(2, 5, 6, 12)

def decompose(x):
    t = x.mean((0, 2), keepdims=True)
    s = x.mean((1, 2), keepdims=True)
    cell = x.mean(2, keepdims=True)
    interaction = cell - t - s
    resid = x - cell
    effects = [t, s, interaction, resid]
    ss = np.array([np.square(np.broadcast_to(e, x.shape)).sum() for e in effects])
    if not np.isclose(ss.sum(), np.square(x).sum()):
        raise ArithmeticError('Orthogonal sums of squares did not close.')
    return effects, ss

def permute(x, effects, obs, n, seed):
    """Sensitivity tests with within-strain/within-treatment exchangeability.

    Treatment labels permute within each strain (30 rows); strain labels within
    each treatment (12 rows). Freedman–Lane interaction test permutes residual
    vectors of the additive model within strain and refits the full design.
    Ordinal sample numbers are not assumed to establish cross-condition day blocks.
    """
    rng = np.random.default_rng(seed)
    s_effect, t_effect = effects[1], effects[0]
    additive = s_effect + t_effect
    redres = x - additive
    xt = x.reshape(2, 30, 12)
    xs = x.transpose(1, 0, 2, 3).reshape(5, 12, 12)
    ri = redres.reshape(2, 30, 12)
    exceed = np.zeros(3, dtype=int)
    for _ in range(n):
        a = np.stack([xx[rng.permutation(30)] for xx in xt]).reshape(x.shape)
        b = np.stack([xx[rng.permutation(12)] for xx in xs]).reshape(5, 2, 6, 12)
        r = np.stack([xx[rng.permutation(30)] for xx in ri]).reshape(x.shape) + additive
        c = r.mean(2) - r.mean((0, 2))[None, :, :] - r.mean((1, 2))[:, None, :]
        stats = [12 * np.square(a.mean((0, 2))).sum(),
                 30 * np.square(b.mean((0, 2))).sum(), 6 * np.square(c).sum()]
        exceed += np.asarray(stats) >= obs[:3] - 1e-11
    return (exceed + 1) / (n + 1)

def make_figure(eff, parts, out):
    te, _, _, res = eff
    centers = te[0, :, 0, :]
    _, d, vt = np.linalg.svd(centers, full_matrices=False)
    load = vt[:2].T
    if (centers[3] - centers[0]) @ load[:, 0] < 0:
        load[:, 0] *= -1
    if load[-2, 1] < 0:
        load[:, 1] *= -1
    pc = d**2 / (d**2).sum() * 100
    contribution = 100 * load[:, 0]**2
    scores = (te + res) @ load
    pd.DataFrame({'endpoint': FIELDS, 'PC1_loading': load[:, 0], 'PC2_loading': load[:, 1],
                  'PC1_contribution_pct': contribution}).to_csv(out/'loadings.csv', index=False)
    sr = []
    for si, strain in enumerate(STRAINS):
        for ti, treatment in enumerate(TREATMENTS):
            for ri, pts in enumerate(scores[si, ti]):
                sr.append([strain, treatment, ri + 1, *pts])
    pd.DataFrame(sr, columns=['strain','treatment','replicate','score1','score2']).to_csv(out/'scores.csv',index=False)
    fig, axs = plt.subplots(2, 2, figsize=(8.2, 6.9))
    fig.subplots_adjust(left=.158, right=.97, bottom=.115, top=.867, hspace=.48, wspace=.65)
    plt.rcParams['font.family'] = 'DejaVu Sans'
    for ax, lab, title in zip(axs.flat,'ABCD',['Treatment scores','Treatment loadings','Contributions to treatment PC1','Partition of total variance']):
        ax.set_title(title, loc='left', fontsize=10, fontweight='bold', pad=9)
        ax.text(-.22,1.1,lab,transform=ax.transAxes,fontsize=13,fontweight='bold')
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=8)
    a,b,c,e = axs.flat
    for ti,(t,color) in enumerate(zip(TREATMENTS,COLORS)):
        points=scores[:,ti].reshape(-1,2)
        a.scatter(*points.T,s=14,edgecolor=color,facecolors='none',linewidth=.8)
        a.scatter(*(centers[ti] @ load),marker='X',s=47,color=color,label=t)
    a.axhline(0,color='#bbb',lw=.6,zorder=0)
    a.set_xlabel(f'PC1 ({pc[0]:.1f}%)', fontsize=9);a.set_ylabel(f'PC2 ({pc[1]:.1f}%)', fontsize=9)
    handles, names = a.get_legend_handles_labels()
    fig.legend(handles, names, ncol=5, loc='upper center', bbox_to_anchor=(.52,.987), fontsize=7.5, frameon=False, columnspacing=1.1, handletextpad=.25)
    for i, (x,y) in enumerate(load):
        b.scatter(x,y,s=16,color='#346779')
        b.annotate(str(i+1),(x,y),xytext=(3,3 if i%2 else -8),textcoords='offset points',fontsize=7)
    b.margins(.19);b.axhline(0,color='#bbb',lw=.6,zorder=0)
    b.set_xlabel('PC1 loading',fontsize=9);b.set_ylabel('PC2 loading',fontsize=9)
    order=np.argsort(contribution)
    c.barh(np.arange(12),contribution[order],color='#3a7787',height=.66)
    c.set_yticks(np.arange(12),[f'{i+1}  {LABELS[i]}' for i in order],fontsize=7)
    c.set_xlim(0,16.5);c.set_xlabel('Contribution (%)',fontsize=9)
    for j,i in enumerate(order):c.text(contribution[i]+.2,j,f'{contribution[i]:.1f}',va='center',fontsize=6.5)
    values=parts
    e.barh(range(4),values,color=['#227565','#5693a4','#a98256','#8b939c'],height=.53)
    e.set_yticks(range(4),['Treatment','Strain','Interaction','Residual'],fontsize=8);e.invert_yaxis()
    e.set_xlim(0,65);e.set_xlabel('Total standardized variance (%)',fontsize=9)
    for j,val in enumerate(values):e.text(val+1,j,f'{val:.2f}',va='center',fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'Figure_7.{ext}',dpi=600,facecolor='white')
    plt.close(fig)
    return pc, contribution

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('asca_input.csv'))
    parser.add_argument('--out',type=Path,default=Path('results'))
    parser.add_argument('--permutations',type=int,default=10000)
    parser.add_argument('--seed',type=int,default=271828)
    args=parser.parse_args()
    if args.permutations<1:parser.error('permutations must be at least 1')
    args.out.mkdir(parents=True,exist_ok=True)
    x=data_cube(args.input);effects,ss=decompose(x)
    p=permute(x,effects,ss,args.permutations,args.seed);parts=100*ss/ss.sum()
    pc,_=make_figure(effects,parts,args.out)
    res={'partition_percent':dict(zip(['treatment','strain','interaction','residual'],parts)),
         'treatment_PC_percent':pc.tolist(),'sensitivity_P':dict(zip(['treatment','strain','interaction'],p)),
         'permutations':args.permutations,'seed':args.seed}
    (args.out/'summary.json').write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps(res,indent=2))
if __name__=='__main__':main()
