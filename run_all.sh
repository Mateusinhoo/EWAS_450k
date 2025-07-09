#!/bin/bash
set -e

echo "Running EWAS for BLOOD samples..."
snakemake --cores 4 --use-conda --latency-wait 30 --configfile config_blood.yml

echo "Running EWAS for BRAIN samples..."
snakemake --cores 4 --use-conda --latency-wait 30 --configfile config_brain.yml

echo "All EWAS runs completed!"
