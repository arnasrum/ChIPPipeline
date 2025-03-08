from abc import ABC, abstractmethod


class PeakCaller(ABC):
    @abstractmethod
    def call_peaks(self,
                   treatment_files: list[str],
                   peak_type: str,
                   paired_end: bool,
                   outfile_name: str,
                   control_files: list[str] = None,
                   args: str = None) -> str:
        pass