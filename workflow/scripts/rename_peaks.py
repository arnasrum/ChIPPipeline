import json
import argparse

def __flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict


def get_macs_input() -> dict[str: dict]:
    with open("config/samples.json") as file:
        samples = json.load(file)
    macs_input = {}
    control_files: list[(str, int, str)] = []
    for entry in __flatten_dict(samples).values():
        replicate = str(entry["replicate"])
        if entry["type"] == "treatment":
            file_name = f"{entry['mark']}_{entry['sample']}"
            if not file_name in macs_input: macs_input[file_name] = {}
            if not replicate in macs_input[file_name]:
                macs_input[file_name][replicate] = {"treatment": [entry['cleanFileName']], "control": [], "peak_type": entry["peak_type"]}
            else:
                macs_input[file_name][replicate]["treatment"].append(entry['cleanFileName'])
        elif entry["type"] == "control":
            control_files.append((entry['sample'], replicate, entry['cleanFileName']))
        else:
            raise Exception(f"Entry type unrecognized for {entry['cleanFileName']}")
    for control in control_files:
        for file_name in filter(lambda file_name: control[0] in file_name, macs_input.keys()):
            macs_input[file_name][control[1]]["control"].append(control[2])

    invalid_input: list[(str, str)] = []
    for file_name in macs_input:
        for missing_pair_replicate in filter(
                lambda replicate: len(macs_input[file_name][replicate]["treatment"]) == 0 or len(
                        macs_input[file_name][replicate]["control"]) == 0, macs_input[file_name].keys()):
            invalid_input.append((file_name, missing_pair_replicate))
    for file, replicate in invalid_input:
        macs_input[file].pop(replicate)
    return macs_input

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