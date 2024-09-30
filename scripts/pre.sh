#!/bin/bash

# Pre-run script, runs after payu sets up the work directory and before the model is run

# Sets the appropriate land use for the model start date.
# The model doesn't vary the land use itself, this is instead done as a
# post-processing step - the old land use values are moved to a new STASH code,
# and new land use values are read in from an external file.

source  /etc/profile.d/modules.sh
module use /g/data/hh5/public/modules
module load conda/analysis3

set -eu

# AFP removed landuse part for piControl run
# Keep pre.sh for adding freshwater 
# 211015
