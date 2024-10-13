#!/bin/bash

# Define the path to the FES_from_Reweighting.py script
FES_SCRIPT="../FES_from_Reweighting.py"


COLVAR_FILE="colvar_opes"
TEMP=300
SIGMA=0.06 # maybe: 3rd column in KERNELS file
NUM_BLOCKS=1 # computes the errors

# Run the FES_from_Reweighting.py script with the specified parameters
python3 $FES_SCRIPT --colvar $COLVAR_FILE --temp $TEMP --sigma $SIGMA --blocks $NUM_BLOCKS
