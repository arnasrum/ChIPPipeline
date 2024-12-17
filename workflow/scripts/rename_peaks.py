import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input

def rename_peaks() -> None:
    macs_input = get_macs_input()
    peak_files = {}
    for sample in macs_input:
        for replicate in macs_input[sample]:
            peak = "narrowPeak" if macs_input[sample][replicate]["peak_type"] == "narrow" else "broadPeak"
            peak_files[f"{sample}_rep{replicate}"] = f"results/macs3/{sample}_rep{replicate}_peaks.{peak}"
    for file_name, path in peak_files.items():
        with open(path, "r") as file:
            lines = file.readlines()
        new_lines = []
        for i, line in enumerate(lines):
            line = line.replace(f"{file_name}_", "")
            new_lines.append(line)
        with open(path, "w") as file:
            file.writelines(new_lines)

if __name__ == "__main__":
    rename_peaks()