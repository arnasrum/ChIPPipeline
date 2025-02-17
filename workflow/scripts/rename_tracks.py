import argparse

def update_titles(tracks_file, new_title_map):
    with open(tracks_file) as file:
        lines = file.readlines()
    replace_lines = {}
    for i, line in enumerate(lines):
        if line == "\n" or line.startswith("#"):
            continue
        if line.startswith("title ="):
            file_type = ""
            line_iter = iter(lines[i:])
            y = 0
            while not file_type.startswith("file_type"):
                file_type = next(line_iter).rstrip("\n")
                y += 1
            file_type = file_type.lstrip("file_type = ")
            replace_lines[i] = f"title = {new_title_map[file_type]}\n"
            if file_type == "bed":
                replace_lines[i + y] = "file_type = narrow_peak\n"
    for num in replace_lines:
        lines[num] = replace_lines[num]
    with open(tracks_file, "w") as file:
        file.writelines(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input file")
    args = parser.parse_args()
    new_title_map = {"bed": "bed", "bigwig": "bigwig"}
    update_titles(args.infile, new_title_map)
