import argparse

def rename_peaks(bed_file: str) -> None:
    """
    Renames peaks within a BED file by removing a specific prefix from each peak name.

    This function reads a BED file, identifies a prefix to remove based on the
    filename (removing path, extension, and "_peaks" suffix), and then replaces
    all occurrences of this prefix within each line. The modified content
    is written back to the original file, overwriting its previous content.

    Args:
        bed_file: The path to the BED file to be modified.
    """
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