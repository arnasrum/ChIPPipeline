from workflow.scripts.trimmers.trimmer import Trimmer
import os

class Cutadapt(Trimmer):
    def trim(self, output_path: str, read1: str, read2:str =None, threads:int=1, args:str=None) -> str:

        file_name = Trimmer.get_file_name(read1, read2)
        read_extensions = [".fastq"]
        if read2:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]
        command = "cutadapt"
        command += f" --cores {threads}"
        if args:
            command += f" {args}"
        command += f" -o {output_path}/{file_name}{read_extensions[0]}"
        if read2:
            command += f" -p {output_path}/{file_name}{read_extensions[1]}"
        command += f" {read1}"
        if read2:
            command += f" {read2}"
        return command