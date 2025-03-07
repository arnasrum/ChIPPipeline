from workflow.scripts.trimmers.trimmer import Trimmer
import os

class Cutadapt(Trimmer):
    def trim(self, output_path: str, read1: str, read2:str =None, threads:int=1, args:str=None) -> str:
        if not os.path.isfile(read1) and False:
            raise FileNotFoundError(f"Provided file; {read1}, does not exist")
        if read2 is not None and not os.path.isfile(read2) and False:
            raise FileNotFoundError(f"Provided file; {read2}, does not exist")
        file_name = os.path.basename(read1.split("/")[-1])
        read_extensions = [""]
        if read2 is not None:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]
        file_name = file_name.strip("_1" "_2")
        command = "cutadapt"
        command += f" --cores {threads}"
        if args is not None:
            command += f" {args}"
        command += f" -o {output_path}/{file_name}{read_extensions[0]}"
        if read2 is not None:
            command += f" -p {output_path}/{file_name}{read_extensions[1]}"
        command += f" {read1}"
        if read2 is not None:
            command += f" {read2}"
        return command