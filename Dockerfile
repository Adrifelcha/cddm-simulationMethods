FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    bash \
    openssh-server

# Create SSH directory and set root password
RUN mkdir /var/run/sshd && \
    echo 'root:adri93' | chpasswd && \
    sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# Install required R packages
RUN R -e "install.packages(c(\
    'scatterplot3d', \
    'plot3D' \
    ), repos='https://cran.rstudio.com/')"

# Create a working directory
WORKDIR /app

# Copy your R code into the container
COPY code/cddm /app/code/cddm

EXPOSE 22

CMD ["/usr/sbin/sshd", "-D"]