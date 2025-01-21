# Simulation Methods for the Circular Drift-Diffusion Model (CDDM)

This repository contains code and materials for comparing different simulation methods for generating data from the Circular Drift-Diffusion Model (CDDM). In this repository, all simulation algorithms are implemented in R.

## Overview

The Circular Drift-Diffusion Model (CDDM) is a sequential sampling model designed for decisions with circular alternatives. This project evaluates four distinct computational methods for generating bivariate data from the CDDM:

1. Direct Emulation: Approximates the underlying stochastic process through discrete time steps
2. Rejection Sampling: Operates in a bounded bivariate space
3. Metropolis-Hastings Algorithm: Uses a bivariate normal proposal distribution
4. Numeric Approximation: Applies probability integral transform using trapezoid integration

Each method generates paired angular choices (in radians) and response times (in seconds), and is evaluated for both computational efficiency and statistical accuracy.

## Repository Structure

.