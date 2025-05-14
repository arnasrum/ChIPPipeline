import argparse

def make_tracks(tracks_file: str, bed: str, bigwig: str, options: dict[str, dict[str, str]]):
    sections = {
                "bed": {
                    "file": bed,
                    "title": "Peaks",
                    "height": "2",
                    "color": "red",
                    "show_labels": "false",
                    "font_size": 10
                ,},
                "bigwig": {
                    "file": bigwig,
                    "title": "Bam Coverage",
                    "height": "2",
                    "color": "#666666",
                    "min_value": "0",
                    "max_value": "5",
                    "number_of_bins": "700",
                    "nans_to_zeros": "true",
                    "summary_method": "mean",
                    "show_data_range": "true",
                    "file_type": "bigwig"
                }
    }
    for section in sections:
        sections[section].update(options[section])
    spacer_height = 0.5
    lines = "[x-axis]"
    for section in sections:
        lines += f"\n[spacer]\nheight = {spacer_height}"
        lines += f"\n[{section}]"
        lines += "\n" + "\n".join(f"{option} = {value}" for option, value in sections[section].items())
    with open(tracks_file, "w") as file:
        file.write(lines)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help="Bed input file")
    parser.add_argument("-b", help="Bigwig input file")
    parser.add_argument("-o", help="Trackfile path file")
    parser.add_argument("-x", help="Max value of bigwig")
    args = parser.parse_args()
    options = {
        "bed": {"title": "Peaks"},
        "bigwig": {"title": "Bam Coverage"},
    }
    make_tracks(args.o, args.a, args.b, options)
