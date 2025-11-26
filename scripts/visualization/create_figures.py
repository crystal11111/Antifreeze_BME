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

    def _collect_type_row(tname, payload):
        """Return (row_dict or None) from a type payload."""
        if not isinstance(payload, dict) or 'summary' not in payload:
            return None
        mm = payload.get('summary', {}).get('mean_metrics', {}) or {}
        sd = payload.get('summary', {}).get('std_dev_metrics', {}) or {}
        # Map fields present in your JSON
        row = {
            'Type': tname.upper(),
            'AUC-ROC': mm.get('auc_roc', np.nan),
            'Accuracy': mm.get('accuracy', np.nan),
            'Precision': mm.get('precision', np.nan),
            'Recall': mm.get('recall_sensitivity', np.nan),
            'Specificity': mm.get('specificity', np.nan),
            'F1': mm.get('f1_score', np.nan),
            'MCC': mm.get('mcc', np.nan),
            # (optional) attach std for the table if you like
            'AUC-ROC (std)': sd.get('auc_roc', np.nan),
            'Accuracy (std)': sd.get('accuracy', np.nan),
            'F1 (std)': sd.get('f1_score', np.nan),
        }
        return row

    if not os.path.exists(results_file):
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.text(0.5, 0.5, 'Run classification.ipynb to generate results',
                ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.savefig(os.path.join(FIGDIR, 'figure3_performance.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(FIGDIR, 'figure3_performance.svg'), bbox_inches='tight')
        print("✓ Figure 3 saved (placeholder)")
        return

    with open(results_file, 'r') as f:
        results = json.load(f)

    # Allowed types (order) — include AFGP if it has metrics
    order = ['type1', 'type2', 'type3', 'type4', 'afgp']
    rows = []
    for t in order:
        if t not in results:
            continue
        payload = results[t]
        # Skip explicit error payloads (like AFGP in your file)
        if isinstance(payload, dict) and 'error' in payload:
            continue
        row = _collect_type_row(t, payload)
        if row:
            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.text(0.5, 0.5, 'No valid metrics in JSON', ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.savefig(os.path.join(FIGDIR, 'figure3_performance.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(FIGDIR, 'figure3_performance.svg'), bbox_inches='tight')
        print("✓ Figure 3 saved (no metrics)")
        return

    # --- Build the figure: AUC-ROC, Accuracy, F1, and a table ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # AUC-ROC
    ax = axes[0, 0]
    ax.bar(df['Type'], df['AUC-ROC'])
    ax.set_ylabel('AUC-ROC')
    ax.set_title('AUC-ROC by AFP Type')
    ax.set_ylim(0, 1)
    ax.grid(axis='y', alpha=0.3)

    # Accuracy
    ax = axes[0, 1]
    ax.bar(df['Type'], df['Accuracy'])
    ax.set_ylabel('Accuracy')
    ax.set_title('Accuracy by AFP Type')
    ax.set_ylim(0, 1)
    ax.grid(axis='y', alpha=0.3)

    # F1
    ax = axes[1, 0]
    ax.bar(df['Type'], df['F1'])
    ax.set_ylabel('F1-Score')
    ax.set_title('F1-Score by AFP Type')
    ax.set_ylim(0, 1)
    ax.grid(axis='y', alpha=0.3)

    # Table (include more columns for detail)
    ax = axes[1, 1]
    ax.axis('off')
    table_cols = ['Type', 'AUC-ROC', 'Accuracy', 'Precision', 'Recall', 'Specificity', 'F1', 'MCC']
    shown = df[table_cols].round(3)
    tbl = ax.table(
        cellText=shown.values,
        colLabels=shown.columns,
        cellLoc='center',
        loc='center',
        bbox=[0, 0, 1, 1]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.6)
    ax.set_title('Performance Summary', pad=16)

    plt.tight_layout()
    # Save CSV for the paper repo
    df.to_csv(os.path.join(FIGDIR, 'figure3_metrics_table.csv'), index=False)
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
