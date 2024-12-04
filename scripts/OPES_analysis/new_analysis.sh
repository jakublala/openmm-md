#!/bin/bash

source activate openmm

# python analysis.py --system CD28-G-UW5P2 --date 241118 --project 241010_FoldingUponBinding  --recompute
# python analysis.py --system CD28-G-UW5P4 --date 241118 --project 241010_FoldingUponBinding  --recompute
# python analysis.py --system CD28-G-UW8P4 --date 241118 --project 241010_FoldingUponBinding
# python analysis.py --system CD28-G-MetaD --date 241120 --project 241010_FoldingUponBinding

# python analysis.py --system Z1-B50W --date 241119-2 --project 241109_INFconstruct
# python analysis.py --system Z1-B50W --date 241120-MetaD --project 241109_INFconstruct
# python deltaG.py --system Z1-B50W --date 241119-2 --project 241109_INFconstruct

export RECOMPUTE=True

# DATE=241202-MetaDPace

python analysis.py --system CD28-G --date 241202-Pace50 --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-G --date 241202-Pace100 --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-G --date 241202-Pace200 --project 241010_FoldingUponBinding --recompute $RECOMPUTE

# python analysis.py --system ASYN-A --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system ASYN-G --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-A --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-B --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-G --date 241128-MetaD --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system CD28-P --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system P53-1 --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system P53-2 --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system P53-E --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system SUMO-1A --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE
# python analysis.py --system SUMO-1C --date $DATE --project 241010_FoldingUponBinding --recompute $RECOMPUTE

# python analysis.py --system Z1-B50W --date 241120-MetaD --project 241109_INFconstruct  --recompute $RECOMPUTE
# python analysis.py --system CD28-G-UW5P2 --date 241121 --project 241010_FoldingUponBinding  --recompute $RECOMPUTE
# python analysis.py --system CD28-G-UW5P4 --date 241121 --project 241010_FoldingUponBinding  --recompute $RECOMPUTE
# python analysis.py --system CD28-G-UW8P4 --date 241121 --project 241010_FoldingUponBinding  --recompute $RECOMPUTE

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