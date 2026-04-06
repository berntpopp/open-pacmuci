FROM condaforge/mambaforge:latest AS base

COPY conda/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

RUN mamba install -n open-pacmuci-tools -c bioconda clair3 && mamba clean -afy

COPY . /app
WORKDIR /app

RUN mamba run -n open-pacmuci-tools pip install .

SHELL ["mamba", "run", "-n", "open-pacmuci-tools", "/bin/bash", "-c"]
ENTRYPOINT ["mamba", "run", "-n", "open-pacmuci-tools", "open-pacmuci"]
CMD ["--help"]
