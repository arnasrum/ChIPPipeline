import json

JSON_LOCATION = "samples.json"

def get_macs_input() -> dict[str: dict]:
    with open(JSON_LOCATION) as file:
        samples = json.load(file)
    macs_input = {}
    control_files: list[(str, int, str)] = []
    for key, entry in __flatten_dict(samples).items():
        replicate = str(entry["replicate"])
        if entry["type"] == "treatment":
            file_name = f"{entry['mark']}_{entry['sample']}"
            if not file_name in macs_input: macs_input[file_name] = {}
            if not replicate in macs_input[file_name]:
                macs_input[file_name][replicate] = {"treatment": [entry['file_name']], "control": [],
                                                    "peak_type": entry['peak_type']}
            else:
                macs_input[file_name][replicate]["treatment"].append(entry['file_name'])
        elif entry["type"] == "control":
            control_files.append((entry['sample'], replicate, entry['file_name']))
        else:
            raise Exception(f"Entry type unrecognized for {entry['file_name']}")
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

def __flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict
