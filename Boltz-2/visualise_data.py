#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.lines import Line2D
import numpy as np

def clean_complex_id(complex_str):
    parts = str(complex_str).split('_')
    if len(parts) >= 3 and parts[-1].startswith('s') and parts[-1][1:].isdigit():
        return "_".join(parts[:-1])
    return complex_str

def get_display_name(col_name: str):
    return col_name.replace("boltz_", "").replace("_", " ").capitalize()

def plot_boltz_metric(df, column, output_dir):
    save_path = Path(output_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    display_name = get_display_name(column)
    
    order = df.groupby("base_complex_id")[column].mean().sort_index(ascending=True).index
    
    plt.figure(figsize=(26, 8))
    sns.boxplot(data=df, x='base_complex_id', y=column, order=order, color='red', width=0.7)
    sns.stripplot(data=df, x='base_complex_id', y=column, order=order, color='black', alpha=0.3, size=2)
    
    plt.title(f"Box Plot of {display_name} (All Seeds Combined)")
    plt.xticks(rotation=90, fontsize=8)
    plt.xlabel("Complex (Aggregated Seeds)")
    plt.ylabel(display_name)
    plt.tight_layout()
    plt.savefig(save_path / f"boxplot_{column}.png", dpi=300)
    plt.close()

def plot_ost_metrics(df, output_dir):
    save_path = Path(output_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    metrics = {
        'ost_rmsd': ('RMSD', 1.5, 2.5, 'orange'),
        'ost_l_rmsd': ('L-RMSD', 2.0, 4.0, 'skyblue'),
        'ost_lddt_lp': ('lDDT_LP', 0.7, 0.5, 'lightgreen')
    }
    
    df['is_rank0'] = df['rank'].astype(str) == '0'
    
    for col, (label, high_thresh, low_thresh, color) in metrics.items():
        if col not in df.columns: 
            continue
        valid_df = df[df[col].notna()]
        if valid_df.empty: 
            continue
        
        order = valid_df.groupby("base_complex_id")[col].mean().sort_index(ascending=True).index
        
        plt.figure(figsize=(26, 8))
        sns.boxplot(data=valid_df, x='base_complex_id', y=col, order=order, color=color, width=0.7, showfliers=False)
        
        sns.stripplot(data=valid_df[~valid_df['is_rank0']], x='base_complex_id', y=col, order=order, color='black', alpha=0.3, size=3, jitter=True)
        sns.stripplot(data=valid_df[valid_df['is_rank0']], x='base_complex_id', y=col, order=order, color='red', marker='D', s=6, edgecolor='black', linewidth=0.5, zorder=10)
        
        plt.axhline(y=high_thresh, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
        plt.axhline(y=low_thresh, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        
        custom_lines = [
            Line2D([0], [0], marker='D', color='w', markerfacecolor='red', markersize=8, label='Rank 0 (All Seeds)'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=5, alpha=0.5, label='Other Ranks'),
            Line2D([0], [0], color='green', linestyle='--', label='High Quality'),
            Line2D([0], [0], color='red', linestyle='--', label='Low Quality')
        ]
        
        plt.legend(handles=custom_lines, loc='upper right', fontsize=10)
        plt.title(f"Aggregated Seed Distribution: {label}")
        plt.xticks(rotation=90, fontsize=9)
        plt.grid(axis='y', linestyle=':', alpha=0.3)
        plt.tight_layout()
        plt.savefig(save_path / f"ost_{label.lower().replace('-', '_')}.png", dpi=300)
        plt.close()

def plot_tm_score(df, output_dir="tm_plots"):
    if 'tm_score' not in df.columns: 
        return
    save_path = Path(output_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    valid_df = df[df['tm_score'].notna()].copy()
    valid_df['is_rank0'] = valid_df['rank'].astype(str) == '0'
    
    order = valid_df.groupby("base_complex_id")["tm_score"].median().sort_index(ascending=True).index
    
    plt.figure(figsize=(26, 8))
    sns.boxplot(data=valid_df, x="base_complex_id", y="tm_score", order=order, color="lightblue", width=0.7, showfliers=False)
    
    sns.stripplot(data=valid_df[~valid_df["is_rank0"]], x="base_complex_id", y="tm_score", order=order, color="black", alpha=0.4, size=4, jitter=True)
    sns.stripplot(data=valid_df[valid_df["is_rank0"]], x="base_complex_id", y="tm_score", order=order, color="red", marker="D", s=8, edgecolor="black", linewidth=0.5, zorder=10)
    
    plt.axhline(y=0.5, color="green", linestyle="--", linewidth=1.5, alpha=0.7, label="TM 0.5")
    plt.axhline(y=0.3, color="orange", linestyle="--", linewidth=1.5, alpha=0.7, label="TM 0.3")
    
    custom_lines = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor="red", markersize=10, label="Rank-0 (All Seeds)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=6, alpha=0.6, label="Other Ranks"),
        Line2D([0], [0], color="green", linestyle="--", label="TM-score 0.5"),
        Line2D([0], [0], color="orange", linestyle="--", label="TM-score 0.3")
    ]
    plt.legend(handles=custom_lines, loc="upper right", fontsize=10)
    plt.title("Aggregated TM-score Distribution per Complex")
    plt.xticks(rotation=90, fontsize=9)
    plt.ylabel("TM-score")
    plt.ylim(0, 1.05)
    plt.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path / "tm_scores_boxplot.png", dpi=300)
    plt.close()

def plot_tm_vs_coverage(df, output_dir="tm_plots"):
    if 'tm_score' not in df.columns or 'tm_coverage' not in df.columns: 
        return
    save_path = Path(output_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    valid_df = df[df["tm_score"].notna() & df["tm_coverage"].notna()]
    
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=valid_df, x="tm_coverage", y="tm_score", hue="base_complex_id", alpha=0.6, s=60, palette="tab20")
    
    plt.axhline(y=0.5, color="green", linestyle="--", alpha=0.5)
    plt.axvline(x=90, color="orange", linestyle="--", alpha=0.5)
    
    plt.xlabel("Residues Aligned (%)")
    plt.ylabel("TM-score")
    plt.title("TM-score vs Alignment Coverage (All Seeds)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path / "tm_vs_coverage.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_top_k_performance(df, output_dir="analysis_plots"):
    save_path = Path(output_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    target_metric = 'ost_lddt_lp'
    confidence_metrics = {
        'boltz_confidence_scores': 'Confidence Score',
        'boltz_complex_plddt_scores': 'Complex pLDDT',
        'boltz_ptm_scores': 'pTM',
        'boltz_iptm_scores': 'ipTM'
    }
    
    n_complexes = df['base_complex_id'].nunique()
    k_vals = list(range(1, n_complexes + 1))
    
    plt.figure(figsize=(10, 8))
    
    oracle_mean = df.groupby('base_complex_id')[target_metric].max().mean()
    plt.axhline(y=oracle_mean, color='black', linestyle='--', label='Oracle')
    
    for col, label in confidence_metrics.items():
        if col not in df.columns:
            continue
        
        mean_bests = []
        for k in k_vals:
            complex_means = df.groupby('base_complex_id')[col].mean()
            top_k_complexes = complex_means.nlargest(k).index
            
            best_lddts = df[df['base_complex_id'].isin(top_k_complexes)].groupby('base_complex_id')[target_metric].max()
            mean_bests.append(best_lddts.mean())
        
        plt.plot(k_vals, mean_bests, marker='o', label=label)
    
    plt.title(f"Mean Best lDDT_LP vs. Top K Complexes (N={n_complexes})")
    plt.xlabel("Top K Complexes")
    plt.ylabel("Mean Best lDDT_LP")
    plt.xlim(0, n_complexes + 1)
    plt.grid(True, alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path / "mean_best_lddt_top_k_complexes.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)
    
    df = pd.read_csv(sys.argv[1])
    df['base_complex_id'] = df['complex_id'].apply(clean_complex_id)
    
    boltz_cols = [c for c in df.columns if c.startswith('boltz_')]
    for col in boltz_cols:
        plot_boltz_metric(df, col, "boltz2_plots")
    
    plot_ost_metrics(df, "ost_plots")
    plot_tm_score(df, "tm_plots")
    plot_tm_vs_coverage(df, "tm_plots")
    plot_top_k_performance(df, "analysis_plots")
