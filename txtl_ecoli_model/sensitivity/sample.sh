#!/bin/bash

# Generate sensitivity analysis samples from the command line
# To run: ./shell-example-sample.sh
# Note that "python -m" runs the __main__ module in the SALib.sample package
# The -B flag prevents the interpreter from compiling bytecode (.pyc files)

python -m SALib.sample.saltelli \
     -n 5000 \
	   -p ./Bounds.txt \
	   -o ./Samples.txt \
	   --delimiter=' ' \
	   --precision=8 \

## Options:
# -p, --paramfile: Your parameter range file (3 columns: parameter name, lower bound, upper bound)
#
# -n, --samples: Sample size.
#				 Number of model runs is N(2D + 2) if calculating second-order indices (default)
#        or N(D + 2) otherwise.
#
# -o, --output: File to output your samples into.
#
# --delimiter (optional): Output file delimiter.
#
# --precision (optional): Digits of precision in the output file. Default is 8.
#
# --max-order (optional): Maximum order of indices to calculate. Choose 1 or 2, default is 2.
#								   Choosing 1 will reduce total model runs from N(2D + 2) to N(D + 2)
#								   Must use the same value (either 1 or 2) for both sampling and analys
