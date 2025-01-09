FROM snakemake/snakemake:latest

WORKDIR /pipeline

COPY ./workflow ./workflow
COPY ./config ./config

RUN conda config --set channel_priority strict
RUN conda config --add channels defaults
RUN pip install snakemake-executor-plugin-slurm
RUN mkdir /input

ENTRYPOINT [ "snakemake" , "--sdm", "conda"]
