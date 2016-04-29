#!/bin/bash

# Perform sensitivity analysis from the command line
# To run: ./shell-example-analyze.sh
# Note that "python -m" runs the __main__ module in the SALib.analyze package
# The -B flag prevents the interpreter from compiling bytecode (.pyc files)

python -m SALib.analyze.sobol \
     -p ./Bounds.txt \
     -Y performance.txt.5000 \
     -c 0 \
     --max-order=2 \
		 --parallel \
		 --processors 8 \
     -r 1000

# Options:
# -p, --paramfile: Your parameter range file (3 columns: parameter name, lower bound, upper bound)
#
# -Y, --model-output-file: File of model output values to analyze
#
# -c, --column (optional): Column of model output file to analyze.
#                If the file only has one column, this argument will be ignored.
#
# --delimiter (optional): Model output file delimiter.
#
# --max-order (optional): Maximum order of indices to calculate.
#               This must match the value chosen during sampling.
#
# -r, --resamples (optional): Number of bootstrap resamples used to calculate confidence intervals on indices. Default 1000.
#
# --parallel (optional): Flag to enable parallel execution with multiprocessing
#
# --processors (optional, int): Number of processors to be used with the parallel option
