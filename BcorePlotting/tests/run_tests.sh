#!/bin/bash

###
### Run Tests associated with BcorePlotting
###

### Usage: sh run_tests.sh /path/to/SOURCE_DIR
#    run_tests.sh must be in the same directory as the test#.R scripts
#    SOURCE_DIR is the path to the directory directly outside of BcorePlotting repo

SCRIPT_DIR=./
SOURCE_DIR=$1 

# Test makeMDSplot and related functions
Rscript $SCRIPT_DIR/test1.R $SOURCE_DIR
rm $SCRIPT_DIR/Rplots.pdf

# Test makeHeatmap and related functions
Rscript $SCRIPT_DIR/test2.R $SOURCE_DIR


