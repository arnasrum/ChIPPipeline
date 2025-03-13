from abc import abstractmethod, ABC
import os

class Trimmer(ABC):
    @abstractmethod
    def trim(self, output_path: str, read1: str, read2:str =None, threads:int=1, args:str=None) -> str:
        pass
    @staticmethod
    def get_file_name(self, read1, read2=None) -> str:
        if not os.path.isfile(read1):
            raise FileNotFoundError(f"Provided file; {read1}, does not exist")
        if read2 is not None and not os.path.isfile(read2):
            raise FileNotFoundError(f"Provided file; {read2}, does not exist")
        file_name = os.path.basename(read1.split("/")[-1].split(".")[0])
        file_name = file_name.strip("_1")
        file_name = file_name.strip("_2")
        return str(file_name)
