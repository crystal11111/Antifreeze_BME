import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import os
from matplotlib.patches import FancyBboxPatch
from Bio import AlignIO, Phylo
from collections import Counter

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA = os.path.join(BASE, 'data')
POS_BY_TYPE = os.path.join(DATA, 'positives', 'by_type')
FIGDIR = os.path.join(DATA, 'figure_for_paper')
RESULTS = os.path.join(BASE, 'model_results')
os.makedirs(FIGDIR, exist_ok=True)

def _conservation_entropy(column):
    residues = [r for r in column if r != '-']
    if not residues: return 0.0
    counts = Counter(residues)
    total = len(residues)
    entropy = -sum((c/total)*np.log2(c/total) for c in counts.values())
    max_entropy = np.log2(min(20, len(counts)))
    return 1 - (entropy/max_entropy if max_entropy>0 else 0.0)

def _smooth(x, w= 9):
    if len(x) < w: return x
    k = w//2
    pad = np.pad(x, (k,k), mode='edge')
    return np.convolve(pad, np.ones(w)/w, mode='valid')

def create_alignment_figure():
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    afp_types = ['type1','type2','type3','type4','afgp']
    for idx, afp_type in enumerate(afp_types):
        ax = axes[idx]
        aln_file = os.path.join(POS_BY_TYPE, f"{afp_type}_c90.aln")
        if os.path.exists(aln_file):
            alignment = AlignIO.read(aln_file, "clustal")
            cons = [_conservation_entropy(alignment[:, i]) for i in range(alignment.get_alignment_length())]
            cons_s = _smooth(np.array(cons), w=9)
            pos = np.arange(len(cons_s))
            ax.fill_between(pos, cons_s, alpha=0.5)
            ax.plot(pos, cons_s, linewidth=1)
            ax.set_xlabel('Position'); ax.set_ylabel('Information content (0–1)')
            ax.set_title(f'{afp_type.upper()} Conservation Profile (n={len(alignment)})')
            ax.set_ylim(0, 1); ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, f'{afp_type} alignment not found', ha='center', va='center', transform=ax.transAxes)
    axes[5].axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'figure2_alignments.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIGDIR, 'figure2_alignments.svg'), bbox_inches='tight')
    print("✓ Figure 2 saved")

def create_performance_figure():
    results_file = os.path.join(RESULTS, 'protein_analysis_metrics.json')
    if os.path.exists(results_file):
        with open(results_file, 'r') as f:
            results = json.load(f)
        types = [t for t in ['type1','type2','type3','type4','afgp'] if t in results]
        data = []
        for t in types:
            r = results[t]
            data.append({
                'Type': t.upper(),
                'ROC-AUC': r.get('mean_roc_auc', 0),
                'PR-AUC': r.get('mean_pr_auc', 0),
                'Accuracy': r.get('mean_accuracy', 0),
                'F1': r.get('mean_f1', 0)
            })
        df = pd.DataFrame(data)
        if df.empty:
            fig, ax = plt.subplots(figsize=(12, 10))
            ax.text(0.5,0.5,'No metrics found in JSON', ha='center', va='center'); ax.axis('off')
        else:
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            axes[0,0].bar(df['Type'], df['ROC-AUC']); axes[0,0].set_ylim(0,1); axes[0,0].set_title('ROC-AUC by AFP Type')
            axes[0,1].bar(df['Type'], df['PR-AUC']); axes[0,1].set_ylim(0,1); axes[0,1].set_title('PR-AUC by AFP Type')
            x = np.arange(len(df)); w = 0.35
            axes[1,0].bar(x-w/2, df['Accuracy'], w, label='Accuracy', alpha=0.7)
            axes[1,0].bar(x+w/2, df['F1'], w, label='F1', alpha=0.7)
            axes[1,0].set_ylim(0,1); axes[1,0].set_xticks(x); axes[1,0].set_xticklabels(df['Type']); axes[1,0].legend()
            axes[1,1].axis('off')
            tbl = axes[1,1].table(cellText=df.round(3).values, colLabels=df.columns,
                                  cellLoc='center', loc='center', bbox=[0,0,1,1])
            tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 2)
        plt.tight_layout()
        df.to_csv(os.path.join(FIGDIR, 'figure3_metrics_table.csv'), index=False)
    else:
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.text(0.5, 0.5, 'Run classification.ipynb to generate results', ha='center', va='center', fontsize=14); ax.axis('off')
    plt.savefig(os.path.join(FIGDIR, 'figure3_performance.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIGDIR, 'figure3_performance.svg'), bbox_inches='tight')
    print("✓ Figure 3 saved")

def create_tree_figure():
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    afp_types = ['type1','type2','type3','type4','afgp']
    for idx, t in enumerate(afp_types):
        ax = axes[idx]
        tf = os.path.join(POS_BY_TYPE, f"{t}_c90.nwk")
        if os.path.exists(tf):
            tree = Phylo.read(tf, "newick")
            try:
                Phylo.draw(tree, axes=ax, do_show=False, 
                          label_func=lambda x: x.name.split('|')[1] if '|' in x.name else x.name,
                          branch_labels=None)
            except TypeError:
                from Bio.Phylo._utils import draw as _draw
                _draw(tree, label_func=lambda c: c.name.split('|')[1] if '|' in c.name else c.name, 
                     do_show=False, axes=ax)
            n_seqs = len(tree.get_terminals())
            ax.set_title(f'{t.upper()} UPGMA Tree (n={n_seqs})', fontweight='bold')
            for sp in ('top','right','left','bottom'): ax.spines[sp].set_visible(False)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlabel('Evolutionary distance', fontsize=9)
        else:
            ax.text(0.5,0.5,f'{t} tree not found', ha='center', va='center', transform=ax.transAxes); ax.axis('off')
    axes[5].axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'figure4_trees.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIGDIR, 'figure4_trees.svg'), bbox_inches='tight')
    print("✓ Figure 4 saved")

def main():
    print("Generating figures for AFP classification paper...\n")
    
    create_alignment_figure()
    create_performance_figure()
    create_tree_figure()
    
    print("\nAll figures generated successfully!")
    print("Figures saved to: data/figure_for_paper/")

if __name__ == "__main__":
    main()
