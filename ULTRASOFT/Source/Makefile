#   **********************************************
#   makefile for Ultrasoft Pseudopotential Program
#   **********************************************
#
# ----------------------------------------------------------
# Set options for fortran compiler if necessary
# ----------------------------------------------------------
#
# Specify name of fortran compiler here:
#
#FC = f77
FC = xlf
#
# Specify any needed compiler options here:
#
FFLAGS= 
#
# -----
# Notes
# -----
#
# It is recommended to choose low optimization levels since the
# code runs quickly anyway.
#
# It is recommended to maintain IEEE compatibility of datafiles
# (see file Doc/IEEE).  For example, for DEC/COMPAQ compilers for
# Alphas, use:
# FFLAGS= -convert big_endian
#
# On most Linux systems you can use "FC = g77" successfully and
# without the need to specify compiler or loader options.  Use of
# Absoft or Intel fortran compilers might require an extra loader
# option (e.g., "-lU77") to point to a library where calls to
# routines like "getarg" can be resolved.  Thus, the use of
# "FC = g77" is recommended as the first recourse when "FC = f77"
# does not work.
#
# The program calls library routines 'edate' and 'itime'.
# On most platforms these are automatically available.
# For IBM rs6000 AIX using the xlf compiler, access these using:
# FFLAGS= -qextname
#
# -----------------------------------------------------
# Normal usage:
# -----------------------------------------------------
#
#   make            Compiles programs and makes runatom.x
#   make install    Installs runatom.x in the ../Bin directory
#   make clean      Cleans up
#
# -----------------------------------------------------
# User does not normally modify below here.
# -----------------------------------------------------

SRCS= runatom.f rwsubs.f initsubs.f aesubs.f pbe.f pssubs.f \
	scgsubs.f numsubs.f relsubs.f bessel.f interpol.f trou.f \
	eissubs.f
OBJS= runatom.o rwsubs.o initsubs.o aesubs.o pbe.o pssubs.o \
	scgsubs.o numsubs.o relsubs.o bessel.o interpol.o trou.o \
	eissubs.o

#------------------------------------------------------------

default: runatom.x

runatom.x:  $(OBJS)
	$(FC) $(OBJS) -o runatom.x

# Note: The $(OBJS) files are made using standard 'implicit rules':
# if, eg, a file 'runatom.o' is needed, it is made via
# 	$(FC) $(FFLAGS) -c -o runatom.o runatom.f

semiclean:
	- rm *.o 	

clean:
	- rm *.o *.x tags

install:
	cp runatom.x ../Bin

#------------------------------------------------------------
#
# Support for editor tag files:
#
tags: 
	ctags $(SRCS)

# For emacs, change 'ctags' to 'etags'

# Usage: after typing
#
#   make tags
#
# you can find a subroutine easily by typing, for example,
#
#   vi -t xctype
#
# The vi editor then automatically opens the file containing subroutine
# 'xctype' and positions the cursor at the start of the subroutine.
