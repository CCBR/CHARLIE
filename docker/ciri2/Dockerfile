FROM nciccbr/ccbr_ubuntu_base_20.04:v6

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install CIRI
WORKDIR /opt2
RUN wget -O CIRI.zip https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip/download && \
    unzip CIRI.zip
ENV PATH="/opt2/CIRI2_v2.0.6/:$PATH"

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
