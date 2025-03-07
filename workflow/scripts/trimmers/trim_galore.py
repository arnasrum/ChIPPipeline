from workflow.scripts.trimmers.trimmer import Trimmer
import os

class TrimGalore(Trimmer):
    def trim(self, output_path, read1, read2=None, threads=1, args=None):
        if not os.path.isfile(read1) and False:
            raise FileNotFoundError(f"Provided file; {read1}, does not exist")
        if read2 is not None and not os.path.isfile(read2) and False:
            raise FileNotFoundError(f"Provided file; {read2}, does not exist")
        file_name = os.path.basename(read1.split("/")[-1].split(".")[0])
        file_name = file_name.strip("_1")
        file_name = file_name.strip("_2")
        read_extensions = [".fastq"]
        if read2 is not None:
            file_name = os.path.commonprefix([file_name, os.path.basename(read2.split("/")[-1])])
            read_extensions = ["_1.fastq", "_2.fastq"]
        command = "trim_galore"
        command += f" -j {threads}"
        command += f" -o {output_path}"
        command += f" --basename {file_name}"
        if read2 is not None:
            command += f" --paired"
        if args is not None:
            command += f" {args}"
        output1 = f"{output_path}/{file_name}{read_extensions[0]}"
        if read2 is not None:
            output2 = f"{output_path}/{file_name}{read_extensions[1]}"
        command += f" {read1}"
        if read2 is not None:
            command += f" {read2}"
        command += f" \nmkdir -p {output_path}"
        if read2 is not None:
            command += f" \nmv {output1.replace("_1.fastq", "")}_val_1.fq {output1}"
            command += f" \nmv {output2.replace("_2.fastq", "")}_val_2.fq {output2}"
        else:
            command += f" \nmv {output1}_trimmed.fq {output1}"
        return command