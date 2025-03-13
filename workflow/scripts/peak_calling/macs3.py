from workflow.scripts.peak_calling.peak_caller import PeakCaller


class Macs3(PeakCaller):
    out_dir = None
    temp_dir = None

    def call_peaks(self,
                   treatment_files: list[str],
                   peak_type: str,
                   paired_end: bool,
                   outfile_name: str,
                   control_files: list[str] = None,
                   args: str = None
    ) -> str:
        command = "macs3 callpeak"
        if args:
            command += f" {args}"
        if self.out_dir:
            command += f" --outdir {self.out_dir}"
        if paired_end:
            command += f" -f BAMPE"
        else:
            command += f" -f BAM"
        command += f" -t {" ".join(treatment_files)} "
        command += f" -n {outfile_name}"
        if control_files:
            command += f" -c {" ".join(treatment_files)} "
        if peak_type == "broad":
            command += " --broad"
            command += (f"\nmv {self.out_dir}/{outfile_name}_peaks.broadPeak " +
             f"{self.out_dir}/{outfile_name}.bed")
        if peak_type == "narrow":
            command += (f"\nmv {self.out_dir}/{outfile_name}_peaks.narrowPeak " +
                        f"{self.out_dir}/{outfile_name}.bed")

        return command