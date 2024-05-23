FROM rocker/shiny
#FROM ubuntu:latest

# Install apt-based dependencies
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y curl git g++ zlib1g-dev make bsdmainutils gawk libopenblas-base wget nano libssl-dev unzip libncurses5-dev python3 python3-pip libbz2-dev liblzma-dev ncbi-tools-bin
RUN apt install -y libqt5svg5-dev python-is-python3 tabix

# Install medaka
RUN pip install medaka opencv-python-headless pyabpoa bgzip

# Make additional tools directory
RUN mkdir /tools

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xvf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure --prefix=/tools && \
    make && \
    make install

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 &&\
    tar -xvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/tools && \
    make && \
    make install

# Install chopper
RUN wget https://github.com/wdecoster/chopper/releases/download/v0.7.0/chopper-linux.zip && \
    unzip chopper-linux.zip && \
    chmod +x chopper && \
    cp chopper /tools/bin

# Install minimap2
RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    make && \
    cp minimap2 /tools/bin

# Install flye
RUN cd /tools && \
    git clone https://github.com/fenderglass/Flye && \
    cd Flye && \
    make

# Install Bandage
RUN mkdir tools/Bandage && \
    cd tools/Bandage && \
    wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_dynamic_v0_8_1.zip && \
    unzip Bandage_Ubuntu_dynamic_v0_8_1.zip && \
    cp Bandage dependencies sample_LastGraph /tools/bin

# Install Prokka
RUN apt-get update -y && \
    apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl && \
    git clone https://github.com/tseemann/prokka.git /tools/prokka && \
    /tools/prokka/bin/prokka --setupdb

# Fix libidn11 compatibility
RUN wget http://mirrors.kernel.org/ubuntu/pool/main/libi/libidn/libidn11_1.33-2.2ubuntu2_amd64.deb && \
    apt install ./libidn11_1.33-2.2ubuntu2_amd64.deb && \
    rm libidn11_1.33-2.2ubuntu2_amd64.deb

# Install HMMER via perl
RUN cpan Bio::SearchIO::hmmer

# Install CheckM
RUN pip install numpy matplotlib pysam checkm-genome

# Install CheckM database files
RUN mkdir /home/checkm_db && \
    cd /home/checkm_db && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar -xvf checkm_data_2015_01_16.tar.gz && \
    rm checkm_data_2015_01_16.tar.gz && \
    checkm data setRoot /home/checkm_db

# Set working directory
WORKDIR /tools

# Install Prodigal
RUN wget https://github.com/hyattpd/Prodigal/archive/refs/tags/v2.6.3.tar.gz && \
    tar -xvf v2.6.3.tar.gz && \
    cd Prodigal-2.6.3 && \
    make install INSTALLDIR=/tools/bin && \
    rm -r /tools/Prodigal-2.6.3 /tools/v2.6.3.tar.gz

# Install pplacer + guppy
RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && \   
    unzip pplacer-linux-v1.1.alpha19.zip && \
    rm /tools/pplacer-linux-v1.1.alpha19.zip

# Install amrfinder
RUN mkdir amrfinder &&\
    cd amrfinder &&\
    URL=`curl -s https://api.github.com/repos/ncbi/amr/releases/latest \
    | grep "browser_download_url.*amrfinder_binaries" \
    | cut -d '"' -f 4` &&\
    curl -sOL "$URL" &&\
    filename=`basename $URL` &&\
    tar xvfz $filename &&\
    # Download latest database
    ./amrfinder -u

# Set ENV
ENV PATH="/tools/bin:/tools/Flye/bin:/tools/prokka/bin:/tools/pplacer-Linux-v1.1.alpha19:/tools/amrfinder:$PATH"
ENV QT_QPA_PLATFORM=offscreen 

# Clean up compressed files and source files
WORKDIR /
RUN rm samtools-1.19.2.tar.bz2 bcftools-1.19.tar.bz2 chopper-linux.zip chopper
RUN rm -r samtools-1.19.2 bcftools-1.19 minimap2 /tools/Bandage
RUN rm /tools/amrfinder/amrfinder_binaries_v3.12.8.tar.gz

# Install R packages
RUN R -e "install.packages(c('stringr', 'dplyr', 'ggplot2', 'plotly', 'shinydashboard', 'shinyalert', 'DT', 'htmlwidgets'))"

# Copy app to /srv/shiny-server/
RUN rm -r /srv/shiny-server/*
COPY app.R /srv/shiny-server/

# Make processing directory
RUN mkdir /home/processing

# Make /home/ writeable to all "users"
RUN chmod -R 777 /home/

# Final setup
WORKDIR /home
EXPOSE 3838