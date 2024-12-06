#!/bin/bash

# Define the path where you want to start JupyterLab
# This example assumes you want to start it in the parent directory of where this script is located
WORKDIR=$(dirname "$(pwd)")

# Define the name of your container image
CONTAINER_IMAGE="ecoli_evolution.sif"

# Start JupyterLab using the Apptainer/Singularity image
# Adjust the path to your Apptainer/Singularity executable as necessary
apptainer exec --bind "${WORKDIR}":/opt/notebooks $CONTAINER_IMAGE \
    jupyter lab --notebook-dir=/opt/notebooks --ip=0.0.0.0 --allow-root  --port 8888

