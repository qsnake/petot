#################################################################
# This is the only section that usually needs to be modified
#
# Set parameters:
#   ATOM = chemical symbol of atom
#   CFG  = atomic or ionic configuration used for generation
#   CFT  = atomic or ionic configuration used for testing

ATOM= ti
CFG= d2s2
CFT= d1s2

#################################################################
#
# Standard operation:
# 
#   1.   make          Generate pseudopotential
#   2.                 Optionally, inspect potential: see below
#   3.   make install  Install pseudopotential file in target directory
#   4.   make clean    Clean up
#
# To inspect the potential after it is generated:
#
#   2a.                View the file *ps*out
#   2b.  make test
#   2c.                View the file *test*out
#       
# For other options, see below
#
#################################################################

# This should be root directory of a7.3.3 package:
ROOT= ../../..

PROG= ${ROOT}/Bin/runatom.x
POT=  ${ROOT}/Pot

AEG= ${ATOM}_ae_${CFG}
AET= ${ATOM}_ae_${CFT}
PSP= ${ATOM}_ps
TST= ${ATOM}_test_${CFT}

default: pseudo

#---------------------------------------------------------------
# run all-electron calculations
#---------------------------------------------------------------

ae: ae_ground ae_test

ae_ground: ${AEG}.ae

ae_test: ${AET}.ae

${AEG}.ae: ${AEG}.adat
	${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy

${AET}.ae: ${AET}.adat
	${PROG} ${AET}.adat ${AET}.out ${AET}.ae ${AET}.atwf ${AET}.logd dummy

#---------------------------------------------------------------
# generate pseudopotential
#---------------------------------------------------------------

pseudo: ${PSP}.uspp

${PSP}.uspp: ${PSP}.adat ${AEG}.ae
	${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp

#---------------------------------------------------------------
# test pseudopotential in test configuration
#---------------------------------------------------------------

test: ${TST}.out

${TST}.out: ${TST}.adat ${AET}.ae ${PSP}.uspp
	${PROG} ${TST}.adat ${TST}.out ${AET}.ae ${TST}.atwf ${TST}.logd ${PSP}.uspp

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
	- rm *_test_*.out *_test_*.atwf *_test_*.logd
