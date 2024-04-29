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

# install circRNA_finder
WORKDIR /opt2
ENV VERSION=v1.2
RUN wget https://github.com/orzechoj/circRNA_finder/archive/refs/tags/v${VERSION}.tar.gz -O circRNA.tar.gz && \
    tar -xzvf circRNA.tar.gz && \
    cd circRNA_finder-${VERSION} && \
    python setup.py install
ENV PATH="/opt2/circRNA_finder-${VERSION}/:$PATH"

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
