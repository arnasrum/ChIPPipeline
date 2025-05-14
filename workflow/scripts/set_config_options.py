import re

def set_module_options(config: dict) -> None:
    '''
        Reads config["modules"] and overwrites the default values.
        Optimally specified when running snakemake through CLI 
        with --config modules="--flag value -f value2"
    '''
    if not config["modules"]: return
    pattern = re.compile(r"(?i)\B--[a-z]*\s|\B-[a-z]\s")
    flags = pattern.findall(config["modules"])
    for i in range(len(flags)):
        argument_start = re.search(flags[i], config["modules"]).end()
        argument_end = len(config["modules"]) if i == len(flags) - 1 else re.search(flags[i + 1], config["modules"]).start()
        argument = config["modules"][argument_start: argument_end]
        match flags[i].rstrip():
            case "-t":
                __set_config_option(config, "trimmer", argument)
            case "--trimmer":
                __set_config_option(config, "trimmer", argument)
            case "-a":
                __set_config_option(config, "aligner", argument)
            case "--aligner":
                __set_config_option(config, "aligner", argument)
            case "-p":
                __set_config_option(config, "peak_caller", argument)
            case "--peak-caller":
                __set_config_option(config, "peak_caller", argument)
            case "-d":
                __set_config_option(config, "duplicate_processor", argument)
            case "--duplicate-processor":
                __set_config_option(config, "duplicate_processor", argument)
            case _:
                raise NotImplementedError(f"Modules keyword; {flags[i]}, flag is not supported")

def set_output_paths(config):
    if config["outdir"] != "" and config["outdir"][-1] != "/": config["outdir"] += "/"
    for config_option in ["benchmarks_path", "results_path", "resources_path", "temp_path", "logs_path"]:
        config[config_option] = config["outdir"] + config[config_option]

def __set_config_option(config:dict, option: str, value: str) -> None:
    config[option] = value.rstrip().lstrip()