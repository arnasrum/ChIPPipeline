from abc import abstractmethod, ABC

class Trimmer(ABC):
    @abstractmethod
    def trim(self, output_path: str, read1: str, read2:str =None, threads:int=1, args:str=None) -> str:
        pass
