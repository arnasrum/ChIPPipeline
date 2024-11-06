FROM snakemake/snakemake:latest

WORKDIR /pipeline

COPY . . 

RUN pip3 install pysradb

ENTRYPOINT [ "snakemake" ]