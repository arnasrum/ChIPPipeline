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
    # Load the data
    data = pd.read_csv(file_path, sep='\t', comment='#')

    # Check the first few lines of the data to understand its structure
    print(data.head())

    # Assuming 'Distance to TSS' is a column name
    # Adjust as necessary based on your file's actual column header.
    distance_col = 'Distance to TSS'
    
    if distance_col not in data.columns:
        raise ValueError(f"Column '{distance_col}' not found in the data. Check the header of your file.")

    # Drop rows with missing values in the distance column
    distances = data[distance_col].dropna()

    counts, bin_edges = np.histogram(distances, bins=100)

    # Determine the threshold
    max_count = np.max(counts)
    threshold = max_count * threshold_fraction

    # Determine start and end bins that exceed the threshold
    start_idx = np.argmax(counts > threshold)
    end_idx = len(counts) - np.argmax(counts[::-1] > threshold) - 1

    # Determine new x-axis limits
    x_limits = (bin_edges[start_idx], bin_edges[end_idx + 1])  # end_idx + 1 to include last bin edge

    # Create the histogram plot
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=100, color='blue', alpha=0.7)
    plt.title('Peak Distribution Relative to TSS')
    plt.xlabel('Distance to TSS (bp)')
    plt.ylabel('Number of Peaks')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim(x_limits)

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Plot saved as {output_path}")

# Example usage
# Assuming the output file is named 'annotatePeaksOutput.txt'

#Assuming the output file is named 'annotatePeaksOutput.txt'


def plot_violin(file_path, output_path='violin_plot.png'):
    """
    Plots a violin plot from HOMER annotatePeaks output showing the distribution
    of distances to TSS across different gene types or categories.

    Parameters:
    - file_path: str, path to the annotatePeaks output file.
    - output_path: str, path where the plot image will be saved.
    """
    # Load the data
    data = pd.read_csv(file_path, sep='\t', comment='#')

    # Print column names to find appropriate columns for plotting
    print(data.columns)

    # Adjust the column names based on your file structure
    distance_col = 'Distance to TSS'  # Update with actual name if different
    category_col = 'Gene Type'        # Using 'Gene Type' as a categorical variable

    if distance_col not in data.columns or category_col not in data.columns:
        raise ValueError("Check that the correct column names are used for 'Distance to TSS' and 'Gene Type'")

    # Drop rows with missing values in the necessary columns
    data_clean = data.dropna(subset=[category_col, distance_col])

    # Create the violin plot
    plt.figure(figsize=(12, 6))
    sns.violinplot(x=category_col, y=distance_col, data=data_clean,
                   scale='width', inner='box', palette='muted')

    # Title and labels
    plt.title('TSS Distance Distribution Across Gene Types')
    plt.xlabel('Gene Type')
    plt.xticks(rotation=45)
    plt.ylabel('Distance to TSS')

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Violin plot saved as {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="Path to the annotatePeaks output file.")
    parser.add_argument("output_path", type=str, help="Path where the plot image will be saved.")
    args = parser.parse_args()
    assert os.path.isfile(args.input_file), "Annotated peaks file is not found."
    plot_peak_distribution(args.input_file, args.output_path + "/peak_distribution.png")
    plot_violin(args.input_file, args.output_path + "/violin_plot.png")

#file = '../H3K4me3_8cell_annotate.txt'
#plot_violin(file)
#plot_peak_distribution(file)
