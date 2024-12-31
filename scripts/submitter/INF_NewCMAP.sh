#! /bin/bash
ssh mmm1486@young.rc.ucl.ac.uk << 'EOF'
    cd /home/mmm1486/projects/openmm-md/scripts/241216_NewCMAP_IFN
    
    # Check if the job is running
    if qstat -j "Z1-B50L10W-CMAP" 2>/dev/null; then
        echo "Job Z1-B50L10W-CMAP is still running. Taking no action."
    else
        echo "Job Z1-B50L10W-CMAP is not running. Starting submit.sh..."
        bash submit.sh
    fi
EOF