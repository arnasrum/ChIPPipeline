from abc import ABC, abstractmethod

class Aligner(ABC):

    @abstractmethod
    def align(self, index, read1, read2=None, threads=1, args=None) -> str:
        pass
    @abstractmethod
    def build_index(self, index: str, output_prefix: str, threads=1, args=None) -> str:
        pass
    @abstractmethod
    def get_index_output(self, prefix, genome) -> list[str]:
        pass
    @abstractmethod
    def get_name(self):
        pass