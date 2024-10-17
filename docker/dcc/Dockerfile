FROM nciccbr/ccbr_ubuntu_base_20.04:v6

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install DCC development version
ENV DCC_VERSION=4418a9a12a7f6e883734459f8117eba5d585b59a
WORKDIR /opt2
RUN wget https://github.com/dieterich-lab/DCC/archive/${DCC_VERSION}.zip -O dcc.zip && \
    unzip dcc.zip && \
    cd DCC-${DCC_VERSION} && \
    pip install numpy==1.21.0 && \
    python setup.py install
ENV PYTHONPATH="/usr/local/lib/python3.8/dist-packages/"
RUN which DCC && DCC -h

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
