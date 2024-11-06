

def getFileNames(sample: str, path: str, ext: str, config: dict) -> list[str]:
    samples = [f"{path}/{sample}_1.{ext}"]
    if config["paired_end"]: samples.append(f"resources/reads/{sample}_2.{ext}")
    return samples

