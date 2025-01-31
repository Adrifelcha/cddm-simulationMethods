FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    bash \
    openssh-server \ 
    git

# Create SSH directory and set root password
RUN mkdir /var/run/sshd && \
    echo 'root:zacatecas' | chpasswd && \
    sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# Install required R packages
RUN R -e "install.packages(c(\
    'scatterplot3d', \
    'plot3D', \
    'circular', \
    'grid', \
    'shape', \
    'geostats', \
    'mvtnorm' \
    ), repos='https://cran.rstudio.com/')"

# Clone the repository
RUN git clone git@github.com:Adrifelcha/cddm-simulationMethods.git

EXPOSE 22

CMD ["/usr/sbin/sshd", "-D"]