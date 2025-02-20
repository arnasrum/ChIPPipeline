import argparse
from collections import defaultdict
import re

def update_options(tracks_file: str):
    with open(tracks_file) as file:
        lines = file.read()
    pattern = r"^\[.*\]$"
    option_map = defaultdict(list)
    sections = re.findall(pattern, lines, re.MULTILINE)
    section = ""
    lines_list = lines.split("\n")
    for i, line in enumerate(lines_list):
        if line.startswith("#") or line.strip() == "":
            continue
        if line in sections:
            section = (i, line)
        else:
            if section:
                option_map[section].append(line)
            else:
                raise Exception("Something went wrong parsing the tracks file")

    old_new_map = {"bed": {"title": "Peaks"}, "bigwig": {"title": "Bam Coverage", "max_value": "5"}}
    replace_options(option_map, old_new_map)
    new_lines = "\n[x-axis]"
    for section in option_map:
        new_lines += f"\n{section[1]}"
        for option in option_map[section]:
            new_lines += f"\n{option}"
    with open(tracks_file, "w") as file:
        file.write(new_lines)

def replace_options(option_map: dict[(int, str), list[str]], old_new_map: dict[str, dict[str, str]]) -> None:
    """
    Replaces generated options by make_tracks_file to make more readable figures.
    (Currently, hardcoded)
    """
    for section in option_map:
        options = option_map[section]
        try:
            file_type = next(filter(lambda option: "file_type" in option, options))
            title = next(filter(lambda option: "title" in option, options))
        except StopIteration:
            continue
        if file_type.split(" = ")[1] == "bed":
            label_option = next(filter(lambda option: "label" in option, options))
            options = [option.replace(label_option, f"show_labels = false") for option in options]
            options = [option.replace(title, f"title = {old_new_map['bed']['title']}") for option in options]
            options = [option.replace(file_type, f"file_type = narrow_peak") for option in options]
        if file_type.split(" = ")[1] == "bigwig":
            options = [option.replace(title, f"title = {old_new_map['bigwig']['title']}") for option in options]
            try:
                max_value = next(filter(lambda option: "max_value" in option, options))
                options = [option.replace(max_value, f"max_value = {old_new_map['bigwig']['max_value']}") for option in options]
            except StopIteration:
                options.append(f"max_value = {old_new_map['bigwig']['max_value']}")
        option_map[section] = options

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input file")
    args = parser.parse_args()
    update_options(args.infile)

if __name__ == "__main__":
    main()
