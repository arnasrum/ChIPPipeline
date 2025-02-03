FROM snakemake/snakemake:v8.20.5

RUN conda install -y conda=23.10.0 snakemake=8.20.5 -c conda-forge -c bioconda  \
    && conda clean -a \
    && pip install snakemake-executor-plugin-slurm \
    && pip install snakemake-executor-plugin-cluster-generic

ENTRYPOINT ["snakemake", "--sdm", "conda", "--conda-frontend", "conda"]

