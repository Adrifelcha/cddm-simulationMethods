# Simulation Methods for the Circular Drift Diffusion Model (CDDM)

This repository contains R code for the implementation of different simulation methods for generating data from the Circular Drift-Diffusion Model (CDDM). The performance of each method is evaluated for both computational efficiency and statistical accuracy.

## Overview

File `cddm_fullReport.Rmd` contains a brief description of the simulation methods explored, and the custom functions build for their implementation. The four simulation methods explored are:

1. **Direct Emulation**: Approximates the underlying stochastic process through discrete time steps
2. **Rejection Sampling**: Operates in a bounded bivariate space
3. **Metropolis-Hastings Algorithm**: Uses a bivariate normal proposal distribution
4. **Probability Integral Transform**: Applies probability integral transform using trapezoid integration

Each method generates paired angular choices (in radians) and response times (in seconds). They are evaluated on the basis of:

- Their execution time
- The difference between the mean angular choice generated and the true underlying drift direction used to generate the data.
- The difference between the empirical cumulative distribution function of the generated data and the theoretical cumulative distribution function approximated through numerical integration.

## Repository Structure

- `code/cddm`: R code specific to the CDDM
- `code/ddm/`: R code specific to the DDM
- `code/figures`: R scripts for generating figures
- `code/general_functions`: R functions for general use
- `code/tests`: R scripts for testing the simulation methods
- `references/`: References and citations
- `tests/`: Placeholder for analysis scripts

# Docker Container for Reproducibility

To ensure reproducibility, we provide a Docker container with all necessary R packages and dependencies.

## Getting Started

### Prerequisites
- Install [Docker](https://docs.docker.com/get-docker/) on your system

### Building and Running the Container

1. Clone this repository:

```bash
git clone git@github.com:Adrifelcha/cddm-simulationMethods.git
```

2. Move to the code directory:

```bash
cd cddm-simulation-methods
```

3. Build the Docker image:

```bash
docker build -t cddm-simulation-methods .
```

3. Run the container:

```bash
docker run -p 2222:22 cddm-simulation-methods
```

This will start the container and expose the SSH service on port 2222. 

4. Connect to the container via SSH:

```bash
ssh -p 2222 root@localhost
```

When prompted, write down the root password: adri93

### Working with the Container

Once connected to the container:
- All required R packages are pre-installed
- The code is available in the `/app/code/cddm` directory
- You can run R scripts and interact with the R environment



