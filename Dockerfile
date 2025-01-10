#FROM condaforge/miniforge3:latest
FROM snakemake/snakemake:v8.26.0


#RUN apt-get -y update && apt-get -y install python3 && apt-get -y install gcc
#RUN pip3 install snakemake==8.26.0 #&& pip3 install pandas
RUN conda config --set channel_priority strict && conda config --add channels defaults

WORKDIR /pipeline
COPY ./workflow ./workflow
COPY ./config ./config
COPY ./profiles ./profiles

#RUN pip install snakemake-executor-plugin-slurm
#RUN pip install snakemake-executor-plugin-cluster-generic
RUN mkdir /input && mkdir /cluster

ENTRYPOINT [ "snakemake", "--sdm", "conda" ]
