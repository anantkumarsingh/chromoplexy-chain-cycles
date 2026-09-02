import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.ma.extras import stack
from scipy import stats
import warnings
import os
import seaborn as sns

warnings.filterwarnings('ignore')

# ── CONFIG ───────────────────────────────────────────────────────────────────
CLINICAL_PATH = "../../../data/clinical/patients_clinical_features.csv"
RESULTS_DIR   = "../../../results/clinical"
os.makedirs(RESULTS_DIR, exist_ok=True)
# ── LOAD DATA ────────────────────────────────────────────────────────────────
clinical = pd.read_csv(CLINICAL_PATH, low_memory=False)
print(f"Patients loaded: {len(clinical)}")

clinical['gleason_score'] = pd.to_numeric(clinical['gleason_score'], errors='coerce')
clinical['gleason_primary'] = pd.to_numeric(clinical['gleason_primary'], errors='coerce')
clinical['gleason_secondary'] = pd.to_numeric(clinical['gleason_secondary'], errors='coerce')
clinical['age_at_initial_pathologic_diagnosis'] = pd.to_numeric(clinical['age_at_initial_pathologic_diagnosis'],
                                                                errors='coerce')

# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1.1A — GLEASON SCORE DISTRIBUTION
# ════════════════════════════════════════════════════════════════════════════

def gleason_distribution(df):
    print("\n" + "=" * 60)
    print("ANALYSIS 1.1A — GLEASON SCORE DISTRIBUTION")
    print("=" * 60)

    gleason_counts = df['gleason_score'].value_counts().sort_index() # count frequency of each gleason score
                                                                    # sort each gleason score in increasing order

    # ── Gleason Score to Grade ────────────────────────────────────
    # Based on standard Gleason Score and Grade Group Comparison
    def gs_grade(row):
        g_score = row['gleason_score']
        gp = row['gleason_primary']
        gs = row['gleason_secondary']
        if g_score == 6:
            return 'Low (6)' # Group 1
        elif g_score == 7 and gp == 3:
            return 'Intermediate-Favorable (3+4)' # Group 2
        elif g_score == 7 and gp == 4:
            return 'Intermediate-Unfavorable (4+3)' # Group 3
        elif g_score == 8:
            return 'High (8)' # Group 4
        elif g_score >= 9:
            return 'Very High (9-10)' # Group 5
        else:
            return 'Unknown'

    df['risk_order'] = df.apply(gs_grade, axis=1) # applies gs_group func row-wise and creates another column 'risk_group'
    risk_counts = df['risk_order'].value_counts() # gets frequency of each gleason grade group

    print("\nClinical risk group distribution:")
    risk_order = [
        'Low (6)',
        'Intermediate-Favorable (3+4)',
        'Intermediate-Unfavorable (4+3)',
        'High (8)',
        'Very High (9-10)'
    ]
    for group in risk_order:
        if group in risk_counts:
            pct = risk_counts[group] / len(df) * 100 # (patients in that group / total patients) x 100
            print(f"  {group:<40} {risk_counts[group]:>4} ({pct:.1f}%)")


    # ── VISUALIZING DATA ────────────────────────────────────
    fig = plt.figure(figsize=(18, 14))
    fig.suptitle('Phase 1.1 — Clinical Cohort Characterization\n'
                 'TCGA Prostate Adenocarcinoma (n=500)',
                 fontsize=16, fontweight='bold', y=0.98)

    gridSpec_layout = gridspec.GridSpec(1, 3, figure=fig,
                                  hspace=0.4, wspace=0.35)

    colors = {6: '#2ecc71', 7: '#f39c12', 8: '#e74c3c', 9: '#c0392b', 10: '#7b241c'}

    # ── Plot 1: Gleason score bar chart ──────────────────────────
    ax1 = fig.add_subplot(gridSpec_layout[0, 0]) # first spot in gridSpec
    scores = gleason_counts.index.tolist()
    counts = gleason_counts.values.tolist()
    bar_colors = [colors.get(s, '#95a5a6') for s in scores]
    bars = ax1.bar(scores, counts, color=bar_colors, edgecolor='white',
                   linewidth=1.5, width=0.6)

    # Add count and percentage labels on bars
    for bar, count in zip(bars, counts):
        pct = count / len(df) * 100
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 2,
                 f'{count}\n({pct:.0f}%)',
                 ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax1.set_xlabel('Gleason Score', fontsize=11)
    ax1.set_ylabel('Number of Patients', fontsize=11)
    ax1.set_title('Gleason Score Distribution', fontsize=12, fontweight='bold')
    ax1.set_xticks(scores)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Shows the gleason distribution in bar chart

    # ----- PRIMARY vs SECONDARY HEATMAP --------------
    ax2 = fig.add_subplot(gridSpec_layout[0,1]) # second spot in gridSpec
    pattern_ct = pd.crosstab(df['gleason_primary'],
                             df['gleason_secondary'])
    sns.heatmap(pattern_ct, annot=True, fmt='d', cmap='YlOrRd', ax=ax2)
    ax2.set_xlabel('Secondary Gleason Pattern', fontsize=11)
    ax2.set_ylabel('Primary Gleason Pattern', fontsize=11)
    ax2.set_title('Primary vs Secondary Gleason Pattern', fontsize=12, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    print("Heatmap for Primary vs Secondary Gleason Score Generated")

    # ----- Stage distribution --------------
    ax3 = fig.add_subplot(gridSpec_layout[0, 2])
    stage_order = ['T2a', 'T2b', 'T2c', 'T3a', 'T3b', 'T4']
    stage_gs_ct = pd.crosstab(df['pathologic_T'], df['risk_order'])
    stage_gs_ct = stage_gs_ct.loc[stage_order]
    stage_gs_ct = stage_gs_ct[risk_order]
    risk_colors = [
        '#2ecc71',  # Low (6)                         → green
        '#f1c40f',  # Intermediate-Favorable (3+4)    → yellow
        '#e67e22',  # Intermediate-Unfavorable (4+3)  → orange
        '#e74c3c',  # High (8)                        → red
        '#c0392b'  # Very High (9-10)                → dark red
    ]

    print("\n" + "=" * 60)
    print("ANALYSIS 1.1B — GLEASON SCORE DISTRIBUTION")
    print("=" * 60)

    print(stage_gs_ct)

    stage_gs_ct.plot(kind = 'bar',stacked = True, rot=0, ax = ax3, color = risk_colors)

    ax3.set_xlabel('Risk Order & Pathologic T Stage', fontsize=11)
    ax3.set_ylabel('Number of Patients', fontsize=11)
    ax3.set_title('Pathologic T Stage Distribution', fontsize=12, fontweight='bold')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    print("Gleason Stage distribution generated")

    plt.savefig(RESULTS_DIR + '/gleason_distribution.png', dpi=150, bbox_inches='tight')
    plt.show()


    return df



if __name__ == "__main__":
    clinical = gleason_distribution(clinical)
