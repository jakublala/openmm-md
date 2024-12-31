#! /bin/bash
ssh mmm1486@young.rc.ucl.ac.uk << 'EOF'
    cd /home/mmm1486/projects/openmm-md/scripts/241216_NewCMAP_FuB
    
    # Check if the job is running
    if qstat -j "CD28-G-LongCMAP" 2>/dev/null; then
        echo "Job CD28-G-LongCMAP is still running. Taking no action."
    else
        echo "Job CD28-G-LongCMAP is not running. Starting longCMAP.sh..."
        bash longCMAP.sh
    fi
EOF