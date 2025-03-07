from abc import ABC, abstractmethod

class Aligner(ABC):

    @abstractmethod
    def align(self, index: str, read1: str, read2: str=None, threads: int=1, args: str=None) -> str:
        pass
    @abstractmethod
    def build_index(self, index: str, output_prefix: str, threads=1, args=None) -> str:
        pass
    @abstractmethod
    def get_index_output(self, prefix: str, genome: str) -> list[str]:
        pass
    @abstractmethod
    def get_name(self) -> str:
        pass