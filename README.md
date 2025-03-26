# ChIP-Seq Pipeline

Snakemake reproducible and extensible chromatin immunoprecipitation sequencing data analysis pipeline.

<img alt="pipeline flowchart" src="docs/images/flow.png" width="450" height="550">


# Usage

## Conda

Installation of conda and snakemake is required to run the pipeline. 

Install conda by following the instructions provided by the official conda [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

Then install the environment by running:

`conda evn create --file environment.yml -p <path>`

or: 

`conda evn create --file environment.yml -n pipeline`

Activate the environment and start using the pipeline by running inside the git repository:

`snakemake -c <num_cores> -C <config_options>`

## Apptainer/Singularity

Using Singularity and Apptainer it is possible run the pipeline without installing Snakemake and conda locally.

Run:

`apptainer pull container_name docker://arnasrum/chippipeline`

Snakemake together conda is available for use in the container.

Run the pipeline by running the following commands inside the git repository:

`./container.sif -c <num_cores> -j <num_jobs> -C <config_options>`

or more explicitly

`apptainer exec container.sif snakemake -c <num_cores> -j <num_jobs> -C <config_options>`

N.B. remember to bind the working directory if it is outside the home directory.


### Running on HPC

If conda is available on your HPC create the environment and run:

`snakemake --jobs <num> --profile profiles/slurm`

If conda is not available use containers.

Pull the container from Dockerhub using Apptainer/Singularity.

Because compute nodes can have limited access to the internet, provide the samples and genome as local files.

Install conda packages before submitting jobs to compute nodes by running:

`apptainer exec <container_name> snakemake --use-conda --conda-create-envs-only`

Then you can submit the pipeline to the workload manager, e.g. using Slurm

`srun --nodes <num_nodes>... apptainer run <container_name> ...`

## Specify samples 


<table>
    <th>Column</th>
    <th>Description</th>
    <th>Required</th>
    <th>Valid Values</th>
    <tr>
        <td>Mark</td>
        <td>Identifier for transcription factor or histone mark for treatment file.</td>
        <td>Yes, only for treatment samples. Leave blank for control samples.</td>
        <td>String</td>
    </tr>
    <tr>
        <td>Sample</td>
        <td>The biological sample or condition the sequences were derived from.</td>
        <td>Yes, used for associating samples with corresponding sample origin.</td>
        <td>String</td>
    </tr>
    <tr>
        <td>Type</td>
        <td>Specifies whether the sample is treatment or control.</td> 
        <td>Yes</td>
        <td>treatment/control</td>
    </tr>
    <tr>
        <td>Peak_type</td>
        <td>Must be narrow/broad, decides if the peak caller will treat the sample as narrow or broad peaks.</td>
        <td>Yes, for treatment files. Can be omitted for control files.</td>
        <td>narrow/broad</td>
    </tr>
    <tr>
        <td>Accession</td>
        <td>GEO accession number for publicly available samples, used to download the sample from GEO if a local file path is not specified. If file_path is specified, this column serves as a prefix for the filename.</td>
        <td>Yes; if file_path is not defined.</td>
        <td>String</td>
    </tr>    
    <tr>
        <td>Genome</td>
        <td>The genome that the sample will the aligned to. Samples will also be separated by the genome they are aligned to.</td>
        <td>Yes</td>
        <td>Genome identifier or path to a FASTA or gzipped FASTA file named after the identifier. Example; hg19 or path/to/genome/hg19.fa.gz</td>
    </tr>
    <tr>
        <td>File_path</td>
        <td>Paths to the reads sample reads. If paired_end is set to true, then two paths must be defined in the column by separating them with ";" character. Files without file_path specified will try </td>
        <td>No, if not specified the pipeline will try to download the sample based on the GEO accession.</td>
        <td>If paired_end is false; path/to/read.fq.gz. If paired_end is true: path/to/read1.fq.gz;path/to/read2.fq.gz </td>
    </tr>
        

</table>


## Editing Config 

Refer to the [configuration documentaion](docs/conf.md).

## Module Options

Overwrite the default configuration by running:

`snakemake -c <num> <rule> --config modules="--flag1 <value1> -f <value2>"`

## Flags

**-t** or **--trim** overwrites the trimmer

**-a** or **--align** overwrites the aligner


