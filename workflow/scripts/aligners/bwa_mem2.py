from workflow.scripts.aligners.aligner import Aligner

class BwaMem2(Aligner):
    def align(self, index, read1, read2=None, threads=1, args=None) -> str:
        command = "bwa-mem2 mem"
        if not read1:
            raise Exception("No read1 provided")
        if not index:
            raise Exception("No index provided")
        if args:
            command += " " + args
        command += f" -t {threads}"
        command += " " + index
        command += f" {read1}"
        if read2:
            command += f" {read2}"
        return command

    def build_index(self, index: str, output_prefix: str, threads=1, args=None) -> str:
        args = args if args else ""
        command = f"bwa-mem2 index {args} -p {output_prefix} {index}"
        return command

    def get_index_output(self, prefix, genome) -> list[str]:
        index_files = [f"{prefix.rstrip('/')}/{genome}.{ext}"
             for ext in ["amb", "ann", "pac", "bwt.2bit.64", "0123"]
        ]
        return index_files

    def get_name(self):
        return "bwa-mem2"