import json

def get_macs_input(json_path) -> dict[str: dict]:
    with open(json_path) as file:
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

def get_all_input(config):
    macs_input = get_macs_input(config["json_path"])
    input_files = []
    path = config["results_path"]
    for key, value in macs_input.items():
        for replicate in value:
            input_files.append(f"{path}/deeptools/{key}_rep{replicate}_heatmap.png")
            input_files.append(f"{path}/deeptools/{key}_rep{replicate}_profile.png")
            input_files.append(f"{path}/pyGenomeTracks/{key}_rep{replicate}.png")
            if str(config["paired_end"]).lower() == "true":
                for file in  macs_input[key][replicate]["control"] + macs_input[key][replicate]["treatment"]:
                    input_files.append(f"{path}/fastqc/{config['trimmer']}/{file}_1_fastqc.html", )
                    input_files.append(f"{path}/fastqc/{config['trimmer']}/{file}_2_fastqc.html")
            else:
                for file in  macs_input[key][replicate]["control"] + macs_input[key][replicate]["treatment"]:
                    input_files.append(f"{path}/fastqc/{config['trimmer']}/{file}_fastqc.html")

    input_files += [*map(lambda sample: f"{config['results_path']}/homer/{sample}/homerResults.html",macs_input.keys())]
    return input_files

def symlink_input(json_path: str, file_name: str) -> None:
    with open(json_path) as file:
        samples = json.load(file)["provided"]
    return next((item for item in samples.values() if item["file_name"] == file_name), None)

def __flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict
