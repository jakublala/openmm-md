### Restart

Doing restart of CD28 and A-synuclein, as it seems it hasn't converged at all yet.


To debug:
```
python ../../src/plumed/main.py \
    --filepath "../../data/241010_FoldingUponBinding/input/CD28/CD28_alpha.pdb" \
    --device_index "nan" \
    --device "cuda" \
    --timestep 2 \
    --mdtime 500 \
    --restart_rfile "../../data/241010_FoldingUponBinding/output/241029/CD28/alpha_1/CD28_alpha.state" \
    --output_dir "../../data/241010_FoldingUponBinding/output/241109/CD28/alpha_1"
```