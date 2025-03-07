from workflow.scripts.aligners.aligner import Aligner

class Bowtie2(Aligner):
    def align(self, index, read1, read2=None, threads=1, args=None) -> str:
        command = "bowtie2 -q"
        if not read1:
            raise Exception("No read1 provided")
        if not index:
            raise Exception("No index provided")

        if args:
            command += " " + args
        if read2 is None:
            command += f" -U {read1}"
        elif read2:
            command += f" -1 {read1} -2 {read2}"
        command += f" -x {index}"
        command += f" --threads {threads}"
        return command

    def build_index(self, index: str, output_prefix: str, threads=1, args=None) -> str:
        if args:
            command = f"bowtie2-build {args} --threads {threads} {index} {output_prefix}"
        else:
            command = f"bowtie2-build --threads {threads} {index} {output_prefix}"
        return command

    def get_index_output(self, prefix, genome) -> list[str]:
        index_files = [f"{prefix.rstrip('/')}/{genome}.{ext}"
             for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2"]
        ]
        return index_files

    def get_name(self):
        return "bowtie2"