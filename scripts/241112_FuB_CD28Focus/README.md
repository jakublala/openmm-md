### Focus on CD28-G

Need to focus on CD28-G to figure out how large the box needs to be.
- How much padding to do. How large the box.
- How step the wall should be.
- How large the barrier should be.


To debug:
```
python ../../src/plumed/main.py \
    --filepath "../../data/241010_FoldingUponBinding/input/CD28/CD28_general.pdb" \
    --device_index "nan" \
    --device "cuda" \
    --timestep 2 \
    --mdtime 500 \
    --output_dir "../../data/241010_FoldingUponBinding/output/241112/CD28/alpha_1" \
    --barrier 200 \
    --padding 5 \
    --upper_wall_at 5
```