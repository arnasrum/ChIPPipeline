from workflow.scripts.trimmers.trimmer import Trimmer
import os

class Fastp(Trimmer):
    def trim(self, output_path, read1, read2=None, threads=1, args=None):
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

        command = "fastp"
        command += f" -w {threads}"
        command += f" -i {read1}"
        command += f" -o {output_path}/{file_name}{read_extensions[0]}"
        if read2 is not None:
            command += f" -I {read2}"
            command += f" -O {output_path}/{file_name}{read_extensions[1]}"
        command += f" -j {output_path}/{file_name}.json"
        command += f" -h {output_path}/{file_name}.html"
        if args is not None:
            command += " " + args
        return command
