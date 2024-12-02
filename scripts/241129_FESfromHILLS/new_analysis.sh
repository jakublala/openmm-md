#!/bin/bash

source activate openmm

export RECOMPUTE=False

DATE=241128-MetaD

python analysis.py --system ASYN-A --date $DATE --project 241010_FoldingUponBinding
python analysis.py --system ASYN-G --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system CD28-A --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system CD28-B --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system CD28-G --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system CD28-P --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system P53-1 --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system P53-2 --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system P53-E --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system SUMO-1A --date $DATE --project 241010_FoldingUponBinding 
python analysis.py --system SUMO-1C --date $DATE --project 241010_FoldingUponBinding 

# python analysis.py --system PQ19-B30L4 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system PQ19-B30L7 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system PQ19-B30L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system PQ19-B40L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system PQ19-B40L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system PQ19-B50L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B30L4 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B30L7 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B30L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B40L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B40L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Q7-B50L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B30L4 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B30L7 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B30L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B40L10 --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B40L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE
# python analysis.py --system Z1-B50L10W --date 241122-Explore --project 241109_INFconstruct --recompute $RECOMPUTE