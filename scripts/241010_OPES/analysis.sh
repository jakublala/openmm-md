#!/bin/bash

# Define the path to the FES_from_Reweighting.py script
FES_SCRIPT="scripts/FES_from_Reweighting.py"


COLVAR_FILE="colvar_opes"
TEMP=300


# Run the FES_from_Reweighting.py script with the specified parameters
python3 $FES_SCRIPT --colvar $COLVAR_FILE --temp $TEMP
