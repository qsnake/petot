############################################################
#                                                          #
#  This is the makefile to create the libatom.a library    #
#                                                          #
############################################################

# For Cray T3E
#F77=f90 
#FCOPTS =  -dp -O2 -c  
#FCOPZE =  -dp -O2 -c

# For IBM SP
#F77=xlf
#FCOPTS =  -O2 -c  
#FCOPZE =  -O2 -c

# For HP 
#F77 = f77
#FCOPTS = +O2 -c
#FCOPZE = +e -c

# For Linux
# 
# F77 = fort77
# F77 = g77
F77 = gfortran
FCOPTS = -g -c -W -Wall -Wextra -Wno-unused-parameter -Wno-unused-label -Wno-unused-variable -fPIC -fbounds-check -fcheck-array-temporaries
FCOPZE = -c

# For IBM and DEC
# 
# F77 = f77
# FCOPTS = -O -c
# FCOPZE = -c

# For DEC alpha compiler
# 
# F77 = fort
# FCOPTS = -O -fast -c
# FCOPSP = -O5 -fast  -c
# FCOPZE = -c

#
#  You shouldn't change after this....
#

MAIN  = charge.o datout.o dat_in.o difnrl.o difrel.o dsolv1.o dsolv2.o etotal.o \
        ext.o input.o orban.o potran.o potrv.o potrw.o prdiff.o \
        pseud2.o pseudb.o pseudk.o pseudo.o pseudt.o pseudv.o \
        velect.o vionic.o wtrans.o pcc_exp.o pcc_poly.o pcc_sin.o
       
XC = atomxc.o exchng.o ggaxc.o ldaxc.o pbexc.o pw92c.o pw92xc.o  pzxc.o 
       
MATH = dottvv.o trnsvv.o dmixp.o gamfnd.o gamfn2.o gaussj.o polcoe.o rtbis2.o \
       rtbist.o sbessj.o splift.o splint.o spliq.o zrbact.o zrbac2.o

EIS  = tinvit.o tridib.o
       
TIME = zedate.o zesec.o


libatom: $(MAIN) $(XC) $(MATH) $(EIS) $(TIME)
	ar -rv libatom.a  $(MAIN) $(XC) $(MATH) $(EIS) $(TIME)      

clean:
	rm -f *.o

#
#
#     Main group of subroutines
#
#

charge.o: charge.f
	 $(F77) $(FCOPTS) charge.f

datout.o: datout.f
	 $(F77) $(FCOPTS) datout.f

dat_in.o: dat_in.f
	 $(F77) $(FCOPTS) dat_in.f

difnrl.o: difnrl.f
	 $(F77) $(FCOPTS) difnrl.f

difrel.o: difrel.f
	 $(F77) $(FCOPTS) difrel.f

dsolv1.o: dsolv1.f
	 $(F77) $(FCOPTS) dsolv1.f

dsolv2.o: dsolv2.f
	 $(F77) $(FCOPTS) dsolv2.f

etotal.o: etotal.f 
	 $(F77) $(FCOPTS) etotal.f 

ext.o: ext.f
	 $(F77) $(FCOPTS) ext.f

input.o: input.f
	 $(F77) $(FCOPTS) input.f

orban.o: orban.f
	 $(F77) $(FCOPTS) orban.f

potran.o: potran.f 
	 $(F77) $(FCOPTS) potran.f 

potrv.o: potrv.f
	 $(F77) $(FCOPTS) potrv.f

potrw.o: potrw.f
	 $(F77) $(FCOPTS) potrw.f 

prdiff.o: prdiff.f
	 $(F77) $(FCOPTS) prdiff.f 

pseud2.o: pseud2.f
	 $(F77) $(FCOPTS) pseud2.f

pseudb.o: pseudb.f
	 $(F77) $(FCOPTS) pseudb.f

pseudk.o: pseudk.f
	 $(F77) $(FCOPTS) pseudk.f

pseudo.o: pseudo.f
	 $(F77) $(FCOPTS) pseudo.f

pseudt.o: pseudt.f
	 $(F77) $(FCOPTS) pseudt.f

pseudv.o: pseudv.f
	 $(F77) $(FCOPTS) pseudv.f

velect.o: velect.f
	 $(F77) $(FCOPTS) velect.f

vionic.o: vionic.f
	 $(F77) $(FCOPTS) vionic.f
 
wtrans.o: wtrans.f
	 $(F77) $(FCOPTS) wtrans.f

pcc_exp.o: pcc_exp.f
	 $(F77) $(FCOPTS) pcc_exp.f

pcc_poly.o: pcc_poly.f
	 $(F77) $(FCOPTS) pcc_poly.f

pcc_sin.o: pcc_sin.f
	 $(F77) $(FCOPTS) pcc_sin.f

#
#
#     XC group of subroutines
#
#

atomxc.o: atomxc.f
	 $(F77) $(FCOPTS) atomxc.f
 
exchng.o: exchng.f 
	 $(F77) $(FCOPTS) exchng.f 
 
ggaxc.o: ggaxc.f 
	 $(F77) $(FCOPTS) ggaxc.f 
 
ldaxc.o: ldaxc.f  
	 $(F77) $(FCOPTS) ldaxc.f  
 
pbexc.o: pbexc.f 
	 $(F77) $(FCOPTS) pbexc.f 
 
pw92c.o: pw92c.f 
	 $(F77) $(FCOPTS) pw92c.f 
 
pw92xc.o: pw92xc.f 
	 $(F77) $(FCOPTS) pw92xc.f 

pzxc.o: pzxc.f 
	 $(F77) $(FCOPTS) pzxc.f 
  
#
#
#     Eispak subroutines
#
#

tinvit.o: tinvit.f
	 $(F77) $(FCOPTS) tinvit.f

tridib.o: tridib.f
	 $(F77) $(FCOPTS) tridib.f

#
#
#     Mathematical subroutines
#
#

dottvv.o: dottvv.f
	 $(F77) $(FCOPTS) dottvv.f

trnsvv.o: trnsvv.f
	 $(F77) $(FCOPTS) trnsvv.f

dmixp.o: dmixp.f
	 $(F77) $(FCOPTS) dmixp.f

gamfnd.o: gamfnd.f
	 $(F77) $(FCOPTS) gamfnd.f

gamfn2.o: gamfn2.f
	 $(F77) $(FCOPTS) gamfn2.f

gaussj.o: gaussj.f
	 $(F77) $(FCOPTS) gaussj.f

polcoe.o: polcoe.f
	 $(F77) $(FCOPTS) polcoe.f

rtbis2.o: rtbis2.f
	 $(F77) $(FCOPTS) rtbis2.f

rtbist.o: rtbist.f
	 $(F77) $(FCOPTS) rtbist.f

sbessj.o: sbessj.f
	 $(F77) $(FCOPTS) sbessj.f

splift.o: splift.f
	 $(F77) $(FCOPTS) splift.f

splint.o: splint.f
	 $(F77) $(FCOPTS) splint.f

spliq.o: spliq.f
	 $(F77) $(FCOPTS) spliq.f

zrbact.o: zrbact.f
	 $(F77) $(FCOPTS) zrbact.f

zrbac2.o: zrbac2.f
	 $(F77) $(FCOPTS) zrbac2.f

#
#
#     Timing subroutines
#
#

zedate.o: zedate.f
	 $(F77) $(FCOPZE) zedate.f

zesec.o: zesec.f
	 $(F77) $(FCOPZE) zesec.f

#
#   END
#
