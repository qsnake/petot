#   **********************************************
#   makefile for Ultrasoft Reformatting Program
#   **********************************************
#
# ----------------------------------------------------------
# Set options for f77 compiler if necessary
# ----------------------------------------------------------
#
# If necessary, uncomment or modify one of the 'FOPTS=' below.
# (And of course, if your compiler is not f77, change that
# too.)
#
# It is recommended to maintain IEEE compatibility of datafiles
# (see file Doc/IEEE).
#
# It is recommended to choose low optimization levels since the
# code runs quickly anyway.
#
#  --- DEC/COMPAQ compilers for Alphas ---
#      (option '-convert big_endian' to set IEEE compatibility)
#
# FOPTS= -convert big_endian 
#
#  --- IBM rs6000 AIX using the xlf compiler ---
#      (option '-qextname' for compatibility with edate and itime )
#
# FOPTS= -qextname
#
# -----------------------------------------------------
# Normal usage:
# -----------------------------------------------------
#
#   make            Compiles reform.f to make reform.x
#   make install    Installs reform.x in the ../Bin directory
#   make clean      Cleans up
#
# -----------------------------------------------------
# User does not normally modify below here.
# -----------------------------------------------------
#
FFLAGS=  $(FOPTS)

default: reform.x

reform.x: reform.f
	f77 reform.f -o reform.x

install:
	cp reform.x ../Bin

clean:
	- rm *.x
