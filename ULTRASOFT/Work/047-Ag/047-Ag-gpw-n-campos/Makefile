#################################################################
# This is the only section that usually needs to be modified
#
# Set parameters:
#   ATOM = chemical symbol of atom
#   CFG  = atomic or ionic configuration used for generation
#   CFT  = atomic or ionic configuration used for testing

ATOM= ag
CFG= ref
CFT1= 1
CFT2= 2
CFT3= 3

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
AET1= ${ATOM}_ae_${CFT1}
AET2= ${ATOM}_ae_${CFT2}
AET3= ${ATOM}_ae_${CFT3}
PSP= ${ATOM}_ps
TST1= ${ATOM}_test_${CFT1}
TST2= ${ATOM}_test_${CFT2}
TST3= ${ATOM}_test_${CFT3}

default: pseudo

#---------------------------------------------------------------
# run all-electron calculations
#---------------------------------------------------------------

ae: ae_ground ae_test

ae_ground: ${AEG}.ae

ae_test: ${AET1}.ae ${AET2}.ae ${AET3}.ae

${AEG}.ae: ${AEG}.adat
	${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy

${AET1}.ae: ${AET1}.adat
	${PROG} ${AET1}.adat ${AET1}.out ${AET1}.ae ${AET1}.atwf ${AET1}.logd dummy

${AET2}.ae: ${AET2}.adat
	${PROG} ${AET2}.adat ${AET2}.out ${AET2}.ae ${AET2}.atwf ${AET2}.logd dummy

${AET3}.ae: ${AET3}.adat
	${PROG} ${AET3}.adat ${AET3}.out ${AET3}.ae ${AET3}.atwf ${AET3}.logd dummy

#---------------------------------------------------------------
# generate pseudopotential
#---------------------------------------------------------------

pseudo: ${PSP}.uspp

${PSP}.uspp: ${PSP}.adat ${AEG}.ae
	${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp

#---------------------------------------------------------------
# test pseudopotential in test configuration
#---------------------------------------------------------------

test: ${TST1}.out ${TST2}.out ${TST3}.out

${TST1}.out: ${TST1}.adat ${AET1}.ae ${PSP}.uspp
	${PROG} ${TST1}.adat ${TST1}.out ${AET1}.ae ${TST1}.atwf ${TST1}.logd ${PSP}.uspp

${TST2}.out: ${TST2}.adat ${AET2}.ae ${PSP}.uspp
	${PROG} ${TST2}.adat ${TST2}.out ${AET2}.ae ${TST2}.atwf ${TST2}.logd ${PSP}.uspp

${TST3}.out: ${TST3}.adat ${AET3}.ae ${PSP}.uspp
	${PROG} ${TST3}.adat ${TST3}.out ${AET3}.ae ${TST3}.atwf ${TST3}.logd ${PSP}.uspp

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
