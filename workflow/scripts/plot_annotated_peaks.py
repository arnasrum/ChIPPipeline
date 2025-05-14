from snakemake.script import snakemake
from snakemake import shell
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def plot_peak_distribution(file_path, output_path='peak_distribution.png', threshold_fraction=0.01):
    """
    Plots the distribution of peaks relative to the TSS from HOMER annotatePeaks output.

    Parameters:
    - file_path: str, path to the annotatePeaks output file.
    - output_path: str, path where the plot image will be saved.
    """
    data = pd.read_csv(file_path, sep='\t', comment='#')

    distance_col = 'Distance to TSS'
    
    if distance_col not in data.columns:
        raise ValueError(f"Column '{distance_col}' not found in the data. Check the header of your file.")

    # Drop rows with missing values in the distance column
    distances = data[distance_col].dropna()

    counts, bin_edges = np.histogram(distances, bins=100)

    # Determine the threshold
    max_count = np.max(counts)
    threshold = max_count * threshold_fraction

    start_idx = np.argmax(counts > threshold)
    end_idx = len(counts) - np.argmax(counts[::-1] > threshold) - 1

    x_limits = (bin_edges[start_idx], bin_edges[end_idx + 1])  # end_idx + 1 to include last bin edge

    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=100, color='blue', alpha=0.7)
    plt.title('Peak Distribution Relative to TSS')
    plt.xlabel('Distance to TSS (bp)')
    plt.ylabel('Number of Peaks')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim(x_limits)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    #print(f"Plot saved as {output_path}")

def plot_violin(file_path, output_path='violin_plot.png'):
    """
    Plots a violin plot from HOMER annotatePeaks output showing the distribution
    of distances to TSS across different gene types or categories.

    Parameters:
    - file_path: str, path to the annotatePeaks output file.
    - output_path: str, path where the plot image will be saved.
    """
    data = pd.read_csv(file_path, sep='\t', comment='#')

    distance_col = 'Distance to TSS'
    category_col = 'Gene Type'

    if distance_col not in data.columns or category_col not in data.columns:
        raise ValueError("Check that the correct column names are used for 'Distance to TSS' and 'Gene Type'")

    data_clean = data.dropna(subset=[category_col, distance_col])

    # Create the violin plot
    plt.figure(figsize=(12, 6))
    sns.violinplot(x=category_col, y=distance_col, data=data_clean,
                   density_norm='width', inner='box')

    plt.title('TSS Distance Distribution Across Gene Types')
    plt.xlabel('Gene Type')
    plt.xticks(rotation=45)
    plt.ylabel('Distance to TSS')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    plot_peak_distribution(snakemake.input[0], snakemake.output["distribution_plot"], snakemake.params["threshold_fraction"])
    plot_violin(snakemake.input[0], snakemake.output["gene_distribution_plot"])
