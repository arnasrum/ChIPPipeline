FROM snakemake/snakemake:latest

WORKDIR /pipeline

COPY ./workflow ./workflow
COPY ./config ./config

RUN conda config --set channel_priority strict
RUN conda config --add channels defaults 

ENTRYPOINT [ "snakemake" , "--use-conda"]
