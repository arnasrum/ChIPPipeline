import json
import argparse
import sys
sys.path.append("workflow/scripts")
from input_scripts import get_macs_input

def __flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict


def rename_peaks(bed_file: str) -> None:
    macs_input = get_macs_input()
    peak_files = {}
    with open(bed_file, "r") as file:
        lines = file.readlines()
    new_lines = []
    for i, line in enumerate(lines):
        line = line.replace(f"{bed_file.split("/")[-1].split(".")[0].replace("_peaks", "")}_", "")
        new_lines.append(line)
    with open(bed_file, "w") as file:
        file.writelines(new_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file")
    args = parser.parse_args()
    rename_peaks(args.input)