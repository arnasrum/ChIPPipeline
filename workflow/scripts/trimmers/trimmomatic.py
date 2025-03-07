from workflow.scripts.trimmers.trimmer import Trimmer
import os

class Trimmomatic(Trimmer):

    run_options = None

    def trim(self, output_path, read1, read2=None, threads=1, args=None) -> str:
        if not os.path.isfile(read1) and False:
            raise FileNotFoundError(f"Provided file; {read1}, does not exist")
        if read2 is not None and not os.path.isfile(read2) and False:
            raise FileNotFoundError(f"Provided file; {read2}, does not exist")
        file_name = os.path.basename(read1.split("/")[-1])
        read_extensions = [".fastq"]
        if read2 is not None:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]
        file_name = file_name.strip("_1" "_2")
        command = "trimmomatic"
        if read2 is not None:
            command += " PE"
        else:
            command += " SE"
        command += f" -threads {threads}"
        command += f" -summary {output_path}/{file_name}_summary.txt"
        if args is not None:
            command += f" {args}"
        command += f" {read1}"
        if read2 is not None:
            command += f" {read2}"
        command += f" {output_path}/{file_name}{read_extensions[0]}"
        unpaired1 = f"{output_path}/{file_name}_unpaired{read_extensions[0]}"
        command += f" {unpaired1}"
        if read2 is not None:
            command += f" {output_path}/{file_name}{read_extensions[1]}"
            unpaired2 = f"{output_path}/{file_name}_unpaired{read_extensions[1]}"
            command += f" {unpaired2}"
        command += f" {self.run_options}"
        command += f"\nrm {unpaired1}"
        if read2 is not None:
            command += f"\nrm {unpaired2}"
        return command