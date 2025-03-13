from workflow.scripts.trimmers.trimmer import Trimmer
import os

class Fastp(Trimmer):
    def trim(self, output_path, read1, read2=None, threads=1, args=None):

        file_name = Trimmer.get_file_name(read1, read2)
        read_extensions = [""]
        if read2:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]

        command = "fastp"
        command += f" -w {threads}"
        command += f" -i {read1}"
        command += f" -o {output_path}/{file_name}{read_extensions[0]}"
        if read2:
            command += f" -I {read2}"
            command += f" -O {output_path}/{file_name}{read_extensions[1]}"
        command += f" -j {output_path}/{file_name}.json"
        command += f" -h {output_path}/{file_name}.html"
        if args:
            command += " " + args
        return command
