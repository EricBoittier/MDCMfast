#!/bin/bash

# Get the MDCMFAST install directory
MDCMFAST_DIR="/home/boittier/Documents/github/MDCMfast"



# Create the directory structure for results
for X in 1-mtp-fit 2-fit-atoms 3-fit-molecule 4-analysis charmm; do
    echo "Creating directory: $MDCMFAST_DIR/results/h2o/$X"
    mkdir "$MDCMFAST_DIR/results/h2o/$X"
done

# Create reference directory if it doesn't exist
mkdir -p "$MDCMFAST_DIR/examples/ref"

echo -e "\nMake sure to place your cube files in:"
echo "$MDCMFAST_DIR/examples/ref/h2o-pot.cube"
echo "$MDCMFAST_DIR/examples/ref/h2o-dens.cube" 
