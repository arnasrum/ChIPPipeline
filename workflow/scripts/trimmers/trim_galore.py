from workflow.scripts.trimmers.trimmer import Trimmer
import os

class TrimGalore(Trimmer):
    def trim(self, output_path, read1, read2=None, threads=1, args=None):

        file_name = Trimmer.get_file_name(read1, read2)
        read_extensions = [".fastq"]
        if read2:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]
        command = "trim_galore"
        command += f" -j {threads}"
        command += f" -o {output_path}"
        command += f" --basename {file_name}"
        if read2:
            command += f" --paired"
        if args:
            command += f" {args}"
        output1 = f"{output_path}/{file_name}{read_extensions[0]}"
        if read2:
            output2 = f"{output_path}/{file_name}{read_extensions[1]}"
        command += f" {read1}"
        if read2:
            command += f" {read2}"
        command += f" \nmkdir -p {output_path}"
        if read2:
            command += f" \nmv {output1.replace("_1.fastq", "")}_val_1.fq {output1}"
            command += f" \nmv {output2.replace("_2.fastq", "")}_val_2.fq {output2}"
        else:
            command += f" \nmv {output1}_trimmed.fq {output1}"
        return command