FROM ubuntu:20.04

# Set environment variables
ENV LANG=C.UTF-8=LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH

# Install necessary system packages
RUN apt-get update --fix-missing && apt-get install -y \
    wget \
    bzip2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_23.10.0-1-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && /opt/conda/bin/conda clean --all -y \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

RUN conda install -y conda=23.10.0 snakemake=8.20.1 -c conda-forge -c bioconda  \
    && conda clean -a 

ENTRYPOINT ["snakemake", "--use-conda", "--conda-frontend", "conda"]

