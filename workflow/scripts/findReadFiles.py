import os
import re
import pathlib
import pandas as pd


def findReads():

    samples = set()
    with open("config/samples.csv", "r") as file:
        for _, row in pd.read_csv(file).iterrows():
            samples = samples.union({sample for sample in row.values[2:]})
    for sample in samples:
        if not re.match(r"GSM[0-9]*", sample):
            path = pathlib.Path(sample)
            files = " ".join(os.listdir(path.parents[0]))
            print(f"sample: {sample}")
            print(f"files: {files}")
            
            print(re.findall(path.stem.split(path.suffix)[0], files))
    

if __name__ == "__main__":
    findReads()