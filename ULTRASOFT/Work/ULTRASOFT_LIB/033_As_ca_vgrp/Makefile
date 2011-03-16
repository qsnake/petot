#################################################################
# This is the only section that usually needs to be modified
#
# Set parameters:
#   ATOM = chemical symbol of atom
#   CFG  = atomic or ionic configuration used for generation

ATOM= as
CFG= s2p3

#################################################################
#
# Standard operation:
# 
#   1.   make          Generate pseudopotential
#   2.                 Optionally, inspect potential by viewing *ps*out
#   3.   make install  Install pseudopotential file in target directory
#   4.   make clean    Clean up
#
# For other options, see below
#
#################################################################

# This should be root directory of a7.3.3 package:
ROOT= ../../..

PROG= ${ROOT}/Bin/runatom.x
POT=  ${ROOT}/Pot

AEG= ${ATOM}_ae_${CFG}
PSP= ${ATOM}_ps

default: pseudo

#---------------------------------------------------------------
# run all-electron calculations
#---------------------------------------------------------------

ae: ae_ground

ae_ground: ${AEG}.ae

${AEG}.ae: ${AEG}.adat
	${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy

#---------------------------------------------------------------
# generate pseudopotential
#---------------------------------------------------------------

pseudo: ${PSP}.uspp

${PSP}.uspp: ${PSP}.adat ${AEG}.ae
	${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp

#---------------------------------------------------------------
# install pseudopotential in target directory
#---------------------------------------------------------------

install:
	cp ${PSP}.uspp ${POT}/`basename \`pwd\``.uspp
	- cp README ${POT}/`basename \`pwd\``.readme
	echo Installed to ${POT}/`basename \`pwd\``.uspp

# The above should have the effect of storing the potential
# as ${POT}/NAME.uspp where NAME is the last segment of the
# pathname of the current directory.

#---------------------------------------------------------------
# clean up
#---------------------------------------------------------------

clean:
	- rm *.out *.ae *.atwf *.logd *.uspp

# to clean up only pseudo outputs, but preserve results of AE runs:
psclean:
	- rm *_ps.out *_ps.atwf *_ps.logd *.uspp
