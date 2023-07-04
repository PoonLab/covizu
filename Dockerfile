FROM ubuntu:latest

RUN apt-get -qy update && \
    apt-get -qqy install python3.10-dev python3.10-distutils binutils build-essential git python3-scipy openmpi-bin openmpi-doc libopenmpi-dev 

ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py

RUN python3.10 /tmp/get-pip.py
RUN /usr/local/bin/pip install biopython mpi4py

ADD http://www.microbesonline.org/fasttree/FastTree.c /root/FastTree.c
RUN gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o /usr/local/bin/fasttree2 /root/FastTree.c -lm && \
    chmod u+x /usr/local/bin/fasttree2

RUN git clone --depth=1 https://github.com/lh3/minimap2.git /root/minimap2 && \
    make -C /root/minimap2 && \
    cp /root/minimap2/minimap2 /usr/local/bin

RUN git clone --depth=1 https://github.com/neherlab/treetime.git /root/treetime && \
    /usr/local/bin/pip install /root/treetime

RUN git clone --depth=1 https://github.com/somme89/rapidNJ.git /root/rapidnj && \
    make -C /root/rapidnj && \
    cp /root/rapidnj/bin/rapidnj /usr/local/bin

RUN git clone https://github.com/PoonLab/covizu.git /root/covizu && \
    git -C /root/covizu submodule init && \
    git -C /root/covizu submodule update