FROM nciccbr/ccbr_ubuntu_base_20.04:v6

RUN apt-get update && apt-get install -y zlib1g zlib1g-dev git build-essential

# Circexplorer2 --> UCSC boostlibraries cufflinks
RUN apt-get install -y libboost1.67-tools-dev

# install conda packages
COPY environment.yml /data2/
ENV CONDA_ENV=clear
RUN mamba env create -n ${CONDA_ENV} -f /data2/environment.yml && \
    echo "conda activate ${CONDA_ENV}" > ~/.bashrc
ENV PATH="/opt2/conda/envs/${CONDA_ENV}/bin:$PATH"
RUN python -m pip install git+https://github.com/YangLab/CLEAR.git
RUN which circ_quant && circ_quant -h

ENV PATH="/opt2:$PATH"
RUN chmod -R a+rX /opt2
WORKDIR /data2
