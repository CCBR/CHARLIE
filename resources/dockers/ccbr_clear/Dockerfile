FROM continuumio/miniconda2

RUN mkdir -p /data2
RUN mkdir -p /opt2

RUN apt-get update && apt-get install -y zlib1g zlib1g-dev git build-essential

# Adding dependencies
# Circexplorer2
# Circexplorer2 --> UCSC dependencies
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/gtfToGenePred
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/genePredToGtf
WORKDIR /opt2
ADD gtfToGenePred /opt2
ADD genePredToGtf /opt2
RUN chmod a+rx gtfToGenePred genePredToGtf
ENV PATH="/opt2:$PATH"
# Circexplorer2 --> bowtie
ADD bowtie-1.1.2.tar.gz /opt2
ENV PATH="/opt2/bowtie-1.1.2:$PATH"
# Circexplorer2 --> UCSC bedtools tophat 
RUN apt-get install -y bedtools
# Circexplorer2 --> UCSC tophat 
ADD tophat-2.1.0.Linux_x86_64.tar.gz /opt2
ENV PATH="/opt2/tophat-2.1.0.Linux_x86_64:$PATH"
# Circexplorer2 --> UCSC boostlibraries cufflinks
RUN apt-get install -y libboost1.67-tools-dev
# wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
ADD cufflinks-2.2.1.Linux_x86_64.tar.gz /opt2
ENV PATH="/opt2/cufflinks-2.2.1.Linux_x86_64:$PATH"
RUN pip install docopt scipy
RUN pip install circexplorer2

# Adding HISAT2
ADD hisat2-2.1.0.tar.gz /opt2
ENV PATH="/opt2/hisat2-2.1.0:$PATH"

# Adding StringTie
ADD stringtie-2.1.4.Linux_x86_64.tar.gz /opt2
ENV PATH="/opt2/stringtie-2.1.4.Linux_x86_64:$PATH"

# INSTALL CLEAR in conda environment
# RUN . /opt/conda/etc/profile.d/conda.sh && \
#     conda create -n clear python=2.7 && \
# 	conda activate clear && \
# 	pip install pybedtools==0.8.1 && \
# 	python -m pip install git+https://github.com/YangLab/CLEAR.git
# RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
#     echo "conda activate clear" >> ~/.bashrc
# # SHELL ["conda", "run", "-n", "clear", "/bin/bash", "-c"]

# # INSTALL CLEAR
WORKDIR /opt2
RUN git clone https://github.com/YangLab/CLEAR
WORKDIR /opt2/CLEAR
RUN python ./setup.py install

RUN chmod -R a+rX /opt2
WORKDIR /data2

MAINTAINER vishal.koparde@nih.gov
