FROM nciccbr/ccbr_ubuntu_base_20.04:v6

RUN apt-get update && apt-get install -y zlib1g zlib1g-dev git build-essential

# Circexplorer2 --> UCSC boostlibraries cufflinks
RUN apt-get install -y libboost1.67-tools-dev

# Adding StringTie
ADD stringtie-2.1.4.Linux_x86_64.tar.gz /opt2
ENV PATH="stringtie-2.1.4.Linux_x86_64:$PATH"

# install conda packages
COPY environment.yml /data2/
ENV CONDA_ENV=clear
RUN mamba env create -n ${CONDA_ENV} -f /data2/environment.yml && \
    echo "conda activate ${CONDA_ENV}" > ~/.bashrc && \
    conda activate ${CONDA_ENV} && \
	python -m pip install git+https://github.com/YangLab/CLEAR.git
ENV PATH="/opt2/conda/envs/${CONDA_ENV}/bin:$PATH"
ENV R_LIBS_USER=/opt2/conda/lib/R/library/


ENV PATH="/opt2:$PATH"
RUN chmod -R a+rX /opt2
WORKDIR /data2
