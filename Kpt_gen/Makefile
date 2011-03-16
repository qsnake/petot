#
# Makefile for stand-alone kpoint generator from xpw
#
SHELL = /bin/sh

PARAM = param.h

FFLAGS = -dp 
PRECISION =
CPP = /opt/ctl/bin/cpp -P -DPARAM_H=\"$(PARAM)\"
FC = f90
COMPRESS = | compress -v
SHAR = shar
TAR = tar cvf -

# Source files:
KPG = kpgen.x
SYM = input.x
CPW = cnvrt.x intpnt.x intsub.x lower.x metric.x ration.x \
	recip.x symchk.x symgen.x warn.x 

SRC = $(KPG) $(SYM) $(CPW)

OBJ=$(SRC:.x=.o)

# Parameter-file dependent objects:
PARAM_DEP = kpgen.o symgen.o intpnt.o intsub.o 

PROGS = kpgen

.SUFFIXES: .x
.x.o:
	$(CPP) $< $(PRECISION) > $*.f
	$(FC) -c $(FFLAGS) $*.f
	rm -f $*.f

kpgen: $(OBJ) $(PARAM)
	$(FC) $(FFLAGS) -o kpgen $(OBJ) 

$(PARAM_DEP): $(PARAM)

list:
	@echo $(SRC)

base: $(KPG) Makefile
	@echo "	Creating compressed tar file..."
	@$(TAR) Makefile $(PARAM) $(KPG) test $(COMPRESS) > kpgbase.tar.Z

tar: $(SRC) Makefile
	@echo "	Creating compressed tar file..."
	@$(TAR) Makefile $(PARAM) $(SRC) test $(COMPRESS) > kpgen.tar.Z

shar: $(SRC) Makefile
	@echo "	Creating compressed shar file..."
	@$(SHAR) Makefile $(PARAM) $(SRC) test test/* $(COMPRESS) > kpgen.sh.Z

collect: $(SRC)
	@(for f in $(SRC); do \
		cat $$f; echo; \
	done;)

clean:
	-rm core *.o *.f a.out
	
realclean: clean rcsclean
	rm $(SYM) $(CPW)

rcsclean:
	-rcsclean $(KPG)

# How to get non-local sources (latest on trunk unless otherwise specified):
$(KPG):
	co $@
$(SYM):
	co $(DSYM)/$@,v
cnvrt.x:
	co $(DCPW)/$@,v
intpnt.x:
	co -r1.2 $(DCPW)/$@,v
intsub.x:
	co $(DCPW)/$@,v
lower.x:
	co $(DCPW)/$@,v
metric.x:
	co $(DCPW)/$@,v
ration.x:
	co $(DCPW)/$@,v
recip.x:
	co $(DCPW)/$@,v
symchk.x:
	co $(DCPW)/$@,v
symgen.x:
	co $(DCPW)/$@,v
warn.x:
	co $(DCPW)/$@,v
