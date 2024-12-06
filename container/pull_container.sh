#!/bin/bash

# Pulls container and replaces old one if it alredy exists
apptainer pull  --force ecoli_evolution.sif  docker://mnytaqw/ecoli_evolution:latest