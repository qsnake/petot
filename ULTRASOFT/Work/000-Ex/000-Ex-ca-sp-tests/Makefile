#################################################################
# This is the only section that usually needs to be modified
#
# Set parameters:
#   ATOM = chemical symbol of atom
#   CFG  = atomic or ionic configuration used for generation
#   CF1  = atomic or ionic configuration used for testing
#   CF2  = another atomic or ionic configuration used for testing

ATOM= ti
CFG= d2s2
CF1= d1s2
CF2= d2s1

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
AE1= ${ATOM}_ae_${CF1}
AE2= ${ATOM}_ae_${CF2}
PSP= ${ATOM}_ps
TS1= ${ATOM}_test_${CF1}
TS2= ${ATOM}_test_${CF2}

default: pseudo

#---------------------------------------------------------------
# run all-electron calculations
#---------------------------------------------------------------

ae: ae_ground ae_test

ae_ground: ${AEG}.ae

ae_test: ${AE1}.ae ${AE2}.ae

${AEG}.ae: ${AEG}.adat
	${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy

${AE1}.ae: ${AE1}.adat
	${PROG} ${AE1}.adat ${AE1}.out ${AE1}.ae ${AE1}.atwf ${AE1}.logd dummy

${AE2}.ae: ${AE2}.adat
	${PROG} ${AE2}.adat ${AE2}.out ${AE2}.ae ${AE2}.atwf ${AE2}.logd dummy

#---------------------------------------------------------------
# generate pseudopotential
#---------------------------------------------------------------

pseudo: ${PSP}.uspp

${PSP}.uspp: ${PSP}.adat ${AEG}.ae
	${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp

#---------------------------------------------------------------
# test pseudopotential in test configuration
#---------------------------------------------------------------

test: ${TS1}.out ${TS2}.out

${TS1}.out: ${TS1}.adat ${AE1}.ae ${PSP}.uspp
	${PROG} ${TS1}.adat ${TS1}.out ${AE1}.ae ${TS1}.atwf ${TS1}.logd ${PSP}.uspp

${TS2}.out: ${TS2}.adat ${AE2}.ae ${PSP}.uspp
	${PROG} ${TS2}.adat ${TS2}.out ${AE2}.ae ${TS2}.atwf ${TS2}.logd ${PSP}.uspp

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
