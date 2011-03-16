# Thereis no default action:

default:

# This makefile is mainly for cleaning up the Bin and Pot directories:

clean-bin:
	- rm Bin/*.x

clean-pot:
	- rm Pot/*.uspp Pot/*.readme

clean:
	- rm Bin/*.x Pot/*.uspp Pot/*.readme

# Also, for convenience, to compile and install the program in one shot:

compile:
	(cd Source; make; make install; make clean)
