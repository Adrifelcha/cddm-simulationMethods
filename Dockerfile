# Use R base image
FROM r-base:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev

# Install required R packages
RUN R -e "install.packages(c('foreach', 'doParallel', 'scatterplot3d', 'plot3D'))"

# Create working directory
WORKDIR /app

# Copy project files
COPY . /app/

# Set default command
CMD ["R"]