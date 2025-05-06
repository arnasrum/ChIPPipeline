# ChIP-Seq Pipeline

Snakemake reproducible and extensible chromatin immunoprecipitation sequencing data analysis pipeline.

<img alt="pipeline flowchart" src="docs/images/flow.png" width="550" height="550">


# Usage

The recommended method to run the pipeline is with Conda. 

Alternatively, it possible to run the pipeline using Apptainer/Singularity containers.


## Conda

Install conda by following the instructions provided by the official conda [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html).

Then install the environment by running:

`conda evn create --file environment.yml`

Activate the environment and start using the pipeline:

`snakemake -c <num_cores> -C <config_options>`

For documentation of config_options, please refer to the [docs](docs/conf.md)

## Apptainer/Singularity

Using Singularity and Apptainer it is possible run the pipeline without installing Snakemake and conda locally.

Run:

`apptainer pull container_name docker://arnasrum/chippipeline`

Snakemake together conda is available for use in the container.

Run the pipeline by running the following commands inside the git repository:

`./container.sif -c <num_cores> -C <config_options>`

or more explicitly:

`apptainer exec container.sif snakemake -c <num_cores> -C <config_options>`

N.B. remember to bind the working directory if it is outside the home directory.


### Running on HPC

If conda is available on your HPC create the environment and run:

`snakemake --jobs <num> --profile profiles/slurm`

If conda is not available, use containers.

Pull the container from Dockerhub using Apptainer/Singularity.

Install conda packages before submitting jobs to compute nodes by running:

`apptainer exec <container_name> snakemake --conda-create-envs-only`

Then you can submit the pipeline to the workload manager, e.g. using Slurm

`srun --nodes <num_nodes>... apptainer run <container_name> ...`

Or edit and submit the provided Slurm [script](slurm.sh):

`sbatch slurm.sh`

## Editing Config 

Refer to the [configuration documentaion](docs/conf.md).

### Preparing Sample Data

Before running the pipeline, the sample sheet needs be prepared, which defines the input data. 

For editing the sample sheet refer to the [sample sheet docs](docs/sample_sheet.md).

