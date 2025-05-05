import argparse
import os
import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, roc_curve, auc
from sklearn.exceptions import UndefinedMetricWarning
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib_venn import venn2
import warnings
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all_svcall", required=True, help="Path to all SV calls file")
    parser.add_argument("--simulated_svcall", required=True, help="Path to simulated SV file")
    parser.add_argument("--output_dir", required=True, help="Directory to save outputs")
    parser.add_argument("--tolerance", type=int, default=500, help="Tolerance for matching positions")
    return parser.parse_args()

def read_sv_file(path):
    return pd.read_csv(path, sep="\t")

def match_variants(df1, df2, tolerance):
    match = []
    for _, row1 in df1.iterrows():
        matched = df2[
            (df2["chrom"] == row1["chrom"]) &
            (df2["svtype"] == row1["svtype"]) &
            (abs(df2["start"] - row1["start"]) <= tolerance) &
            (abs(df2["end"] - row1["end"]) <= tolerance)
        ]
        match.append(1 if not matched.empty else 0)
    return match

def evaluate(tool_df, sim_df, tolerance, total_possible_sites=48000):
    tp_matches = match_variants(tool_df, sim_df, tolerance)
    TP = sum(tp_matches)
    FP = len(tp_matches) - TP

    fn_matches = match_variants(sim_df, tool_df, tolerance)
    FN = fn_matches.count(0)

    P = TP + FN
    N = total_possible_sites - P
    TN = max(N - FP, 0)

    y_pred = [1]*TP + [1]*FP + [0]*FN + [0]*TN
    y_true = [1]*TP + [0]*FP + [1]*FN + [0]*TN

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        accuracy = accuracy_score(y_true, y_pred)
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = auc(fpr, tpr)

    return {
        "TP": TP, "FP": FP, "FN": FN, "TN": TN,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "Accuracy": accuracy,
        "AUC": roc_auc,
        "FPR_list": fpr.tolist(),
        "TPR_list": tpr.tolist()
    }

def plot_roc(tool_metrics, output_dir, tolerance):
    plt.figure()
    for tool, metrics in tool_metrics.items():
        plt.plot(metrics["FPR_list"], metrics["TPR_list"], label=f"{tool} (AUC={metrics['AUC']:.2f})")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve for all tools ({tolerance} bp tolerance)")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"roc_curve_tol{tolerance}.png"))
    plt.close()

def plot_confusion_matrices(metrics, output_dir, tolerance):
    for tool, m in metrics.items():
        cm = np.array([[m["TP"], m["FP"]], [m["FN"], m["TN"]]])
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)

        for i in range(2):
            for j in range(2):
                val = cm[i, j]
                color = plt.cm.Blues(val / cm.max())
                rect = patches.Rectangle((j, 1 - i), 1, 1, facecolor=color, edgecolor='black', linewidth=2)
                ax.add_patch(rect)
                ax.text(j + 0.5, 1 - i + 0.5, str(val),
                        va='center', ha='center', fontsize=14, weight='bold',
                        color='white' if val > cm.max()/2 else 'black')

        ax.text(-0.5, 1.5, "Actual: SV", va='center', ha='left', fontsize=10)
        ax.text(-0.5, 0.5, "Actual: No SV", va='center', ha='left', fontsize=10)
        ax.text(0.5, 2.05, "Predicted: SV", va='bottom', ha='center', fontsize=10)
        ax.text(1.5, 2.05, "Predicted: No SV", va='bottom', ha='center', fontsize=10)

        plt.title(f"Confusion Matrix: {tool} ({tolerance} bp tolerance)", fontsize=12, pad=40)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"confusion_matrix_{tool}_tol{tolerance}.png"), dpi=300)
        plt.close()

def sv_size_distribution(all_df, sim_df, output_dir, tolerance):
    size_bins = [(20, 100), (100, 1000), (1000, 10000), (10000, 100000), (100000, np.inf)]
    bin_labels = ["20-100", "100-1000", "1000-10000", "10000-100000", ">100000"]

    tool_sizes = {}
    for tool in all_df["tool"].unique():
        tool_df = all_df[all_df["tool"] == tool]
        tool_sizes[tool] = {label: 0 for label in bin_labels}

        for (low, high), label in zip(size_bins, bin_labels):
            count = 0
            for _, row1 in tool_df.iterrows():
                if not (low <= row1["length"] < high):
                    continue
                matched = sim_df[
                    (sim_df["chrom"] == row1["chrom"]) &
                    (abs(sim_df["start"] - row1["start"]) <= tolerance) &
                    (abs(sim_df["end"] - row1["end"]) <= tolerance)
                ]
                if not matched.empty:
                    count += 1
            tool_sizes[tool][label] = count

    sns.set(style="whitegrid")
    plt.figure(figsize=(12, 7))
    palette = sns.color_palette("Set2", n_colors=len(tool_sizes))

    for (tool, sizes), color in zip(tool_sizes.items(), palette):
        plt.plot(bin_labels, list(sizes.values()), label=tool, marker='o', linewidth=2, markersize=6, color=color)

    plt.xlabel("SV Size Range (bp)", fontsize=14)
    plt.ylabel("Count of SVs", fontsize=14)
    plt.title(f"SV Size Distribution Across Tools ({tolerance} bp tolerance)", fontsize=16, weight='bold')
    plt.xticks(rotation=45, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title="Tool", fontsize=11, title_fontsize=12, loc="upper right")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"sv_size_distribution_tol{tolerance}.png")
    plt.savefig(output_path, dpi=300)
    plt.close()

def sv_type_distribution(all_df, sim_df, output_dir, tolerance):
    svtypes = ["INS", "DEL"]
    type_counts = {}

    for tool in all_df["tool"].unique():
        tool_df = all_df[all_df["tool"] == tool]
        type_counts[tool] = {}

        for svtype in svtypes:
            sub_df = tool_df[tool_df["svtype"] == svtype]
            count = 0
            for _, row in sub_df.iterrows():
                matched = sim_df[
                    (sim_df["chrom"] == row["chrom"]) &
                    (sim_df["svtype"] == row["svtype"]) &
                    (abs(sim_df["start"] - row["start"]) <= tolerance) &
                    (abs(sim_df["end"] - row["end"]) <= tolerance)
                ]
                if not matched.empty:
                    count += 1
            type_counts[tool][svtype] = count

    df = pd.DataFrame(type_counts).T.fillna(0)

    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 7))
    df.plot(kind="bar", stacked=True, color=sns.color_palette("Set2"), edgecolor='black')

    plt.xlabel("SV Calling Tool", fontsize=14)
    plt.ylabel("Count of Matched SVs", fontsize=14)
    plt.title(f"SV Type Distribution (INS vs DEL, {tolerance} bp tolerance)", fontsize=12, weight='bold')
    plt.xticks(rotation=45, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title="SV Type", fontsize=11, title_fontsize=12, loc="upper right")
    plt.grid(axis='y', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"sv_type_distribution_tol{tolerance}.png")
    plt.savefig(output_path, dpi=300)
    plt.close()

def venn_diagrams(all_df, sim_df, output_dir, tolerance):
    for tool in all_df["tool"].unique():
        tool_df = all_df[all_df["tool"] == tool]
        samples = tool_df["sample"].unique()
        sim_tool_df = sim_df[sim_df["sample"].isin(samples)]

        def get_keys(df):
            return set(
                df["chrom"] + ":" +
                df["start"].astype(str) + "-" +
                df["end"].astype(str) + ":" +
                df["svtype"]
            )

        def get_matched_keys(df1, df2):
            matched_keys = set()
            for _, row in df1.iterrows():
                matched = df2[
                    (df2["chrom"] == row["chrom"]) &
                    (df2["svtype"] == row["svtype"]) &
                    (abs(df2["start"] - row["start"]) <= tolerance) &
                    (abs(df2["end"] - row["end"]) <= tolerance)
                ]
                if not matched.empty:
                    key = f"{row['chrom']}:{row['start']}-{row['end']}:{row['svtype']}"
                    matched_keys.add(key)
            return matched_keys

        detected = get_keys(tool_df)
        simulated = get_keys(sim_tool_df)
        intersection = get_matched_keys(tool_df, sim_tool_df)

        only_detected = detected - intersection
        only_simulated = simulated - intersection

        plt.figure(figsize=(6, 6))
        v = venn2(subsets=(len(only_detected), len(only_simulated), len(intersection)),
                  set_labels=(f"{tool}", "Simulated"))

        v.get_patch_by_id('10').set_color("#66c2a5") 
        v.get_patch_by_id('01').set_color("#fc8d62")  
        v.get_patch_by_id('11').set_color("#8da0cb") 

        for subset in ('10', '01', '11'):
            patch = v.get_patch_by_id(subset)
            if patch:
                patch.set_alpha(0.8)
                patch.set_edgecolor('black')

        for text in v.set_labels:
            text.set_fontsize(12)
            text.set_fontweight('bold')

        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

        plt.title(f"Venn Diagram: {tool} ({tolerance} bp tolerance)", fontsize=14, weight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{tool}_venn_diagram_tol{tolerance}.png"), dpi=300)
        plt.close()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    all_df = read_sv_file(args.all_svcall)
    sim_df = read_sv_file(args.simulated_svcall)

    if "length" not in all_df.columns:
        all_df["length"] = all_df["end"] - all_df["start"]
    if "length" not in sim_df.columns:
        sim_df["length"] = sim_df["end"] - sim_df["start"]

    tools = all_df['tool'].unique()
    metrics = {}

    for tool in tools:
        print(f"Processing: {tool}")
        tool_df = all_df[all_df["tool"] == tool]
        samples = tool_df["sample"].unique()
        sim_tool_df = sim_df[sim_df["sample"].isin(samples)]
        metrics[tool] = evaluate(tool_df, sim_tool_df, args.tolerance)

    metrics_df = pd.DataFrame.from_dict(metrics, orient="index")
    metrics_df.reset_index().rename(columns={'index': 'tool'}).to_csv(
        os.path.join(args.output_dir, "sv_comparison_metrics.tsv"), sep="\t", index=False
    )

    plot_roc(metrics, args.output_dir, args.tolerance)
    sv_size_distribution(all_df, sim_df, args.output_dir, args.tolerance)
    sv_type_distribution(all_df, sim_df, args.output_dir, args.tolerance)
    plot_confusion_matrices(metrics, args.output_dir, args.tolerance)
    venn_diagrams(all_df, sim_df, args.output_dir, args.tolerance)

if __name__ == "__main__":
    main()
