from workflow.scripts.aligners.aligner import Aligner

class STAR(Aligner):

    temp_dir = ""
    output_prefix = ""

    def align(self, index, read1, read2=None, threads=1, args=None) -> str:
        command = "STAR --readFilesType Fastx"
        if not read1:
            raise Exception("No read1 provided")
        if not index:
            raise Exception("No index provided")
        if args:
            command += " " + args
        command += f" --genomeDir {index}"
        command += f" --runThreadN {threads}"
        command += " " + index
        command += f" --readFilesIn {read1}"
        if read2:
            command += f" {read2}"
        command += f" --outFileNamePrefix {self.output_prefix}"
        command += " --outStd SAM"
        return command
    
    def build_index(self, index: str, output_prefix: str, threads=1, args=None) -> str:
        args = args if args else ""
        command = f"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output_prefix} --genomeFastaFiles {index} {args}"
        return command

    def get_index_output(self, prefix, genome) -> list[str]:
        index_files = [f"{prefix.rstrip('/')}/{genome}/{file}"
             for file in ["SA", "SAindex", "Genome"]
        ]
        return index_files

    def get_name(self):
        return "STAR"