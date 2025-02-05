import csv
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse
from os import listdir

def sample_tool_map(paths):
    samples = defaultdict(set)
    for path in paths:
        for file in listdir(path):
            tool = path.split("/")[-1]
            samples[file].add((tool, get_median(path + "/" + file)))
    print(samples)
    return samples

def get_median(path):
    with open(path) as f:
        rows = [row for row in csv.reader(f, delimiter="\t")]
    time_values = [float(item[0]) for item in rows[1:]]
    time_values.sort()
    return calculate_median(time_values)

def calculate_median(values):
    n = len(values)
    mid = n // 2
    if n % 2 == 0:  # even number of items
        median = (values[mid - 1] + values[mid]) / 2.0
    else:  # odd number of items
        median = values[mid]
    return median


def plot(data):
    # Extract data for plotting
    for filename, methods in data.items():
        labels = []
        values = []

        # Collect labels and values
        for method, value in methods:
            labels.append(method)
            values.append(value)

        # Plotting the histogram
        plt.figure()  # Create a new figure for each plot
        plt.bar(labels, values)  # Create a bar chart
        plt.title(f"Histogram for {filename}")
        plt.xlabel('Tool')
        plt.ylabel('TIME')
        plt.xticks(rotation=45)  # Rotate x-axis labels for readability
        plt.tight_layout()  # Adjust layout to fit labels
        # Optionally save each plot as an image file
        # plt.savefig(f"{filename}_histogram.png")
        plt.show()


def make_paths(directory, tools):
    tool_directory = directory + "/benchmarks"
    if tools:
        return [tool_directory + f"/{tool}" for tool in tools[0]]
    else:
        return [tool_directory + f"/{tool}" for tool in listdir(tool_directory)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", action="append", nargs="+", help="The path to the benchmark to parse")
    parser.add_argument("-t", action="append", nargs="+", help="Tools which should be parsed")
    args = parser.parse_args()
    paths = args.directory[0]
    for path in args.directory[0]:
        tool_paths = make_paths(path, args.t)

        samples = sample_tool_map(tool_paths)
        plot(samples)