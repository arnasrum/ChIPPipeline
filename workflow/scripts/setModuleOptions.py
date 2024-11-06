import re

def setModuleOptions(config: dict) -> None:
    '''
        Reads config["modules"] and overwrites the default values.
        Optimally specified when running snakemake through CLI 
        with --config modules="--flag value -f value"
    '''
    # --flag- is not caught as a flag in regex
    # multi character flag with single dash is read as value
    pattern = re.compile(r"(?i)\B--[a-z]*\s|\B-[a-z]\s")
    flags = pattern.findall(config["modules"])
    for i in range(len(flags)):
        argumentStart = re.search(flags[i], config["modules"]).end()
        argumentEnd = len(config["modules"]) if i == len(flags) - 1 else re.search(flags[i + 1], config["modules"]).start()
        argument = config["modules"][argumentStart: argumentEnd]
        match flags[i].rstrip():
            case "-t":
                __setConfigOption(config, "trimmer", argument)
            case "--trim":
                __setConfigOption(config, "trimmer", argument)
            case "-a":
                __setConfigOption(config, "aligner", argument)
            case "--align":
                __setConfigOption(config, "aligner", argument)
            case "-g":
                __setConfigOption(config, "genome", argument)
            case "--genome":
                __setConfigOption(config, "genome", argument)
            case _:
                raise NotImplementedError(f"{flags[i]} flag is not supported")

def __setConfigOption(config:dict, option: str, value: str) -> None:
    config[option] = value.rstrip().lstrip()