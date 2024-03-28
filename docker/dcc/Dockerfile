FROM nciccbr/ccbr_ubuntu_base_20.04:v6

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install conda packages
COPY environment.txt /data2/
RUN mamba install -c conda-forge -c bioconda --file /data2/environment.txt
ENV R_LIBS_USER=/opt2/conda/lib/R/library/

# install DCC
ENV DCC_VERSION=v0.5.0
RUN wget https://github.com/dieterich-lab/DCC/archive/refs/tags/v0.5.0.tar.gz && \
    tar -xzf DCC-${DCC_VERSION}.tar.gz && \
    cd DCC-${DCC_VERSION} && \
    python setup.py install
RUN which DCC && DCC -h

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
