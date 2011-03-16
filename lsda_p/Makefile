#Makefile for Escan

FC=f90
FFLAGS= -dp -O3 -Ounroll2
#FFLAGS= -dp -Ounroll2
#FFLAGS= -dp -g 
LIB= -lapp -lc -lm -lu -lf -lsci -lfi -lsma -lcomm

MODULES=load_data.f \
fft_data.f \
data.f

SRCL0=mainMV.f

SRCL1=moment.f

SRCL=Hpsi_comp.f\
CG_comp.f\
CG_new.f\
getpot2.f\
getpot3.f\
getpot4.f\
getpot4_force.f\
getpot5.f\
getpot5_force.f\
GGAPBE.f\
d3fft_comp.f\
cfft.f\
cfftd.f\
diag_comp.f\
djacobi.f\
d3fft_real2.f\
fftprep_comp.f\
fftprep_real2.f\
fwdcpfft_comp.f\
fwdcpfft2.f\
invcpfft_comp.f\
invcpfft2.f\
gen_G_comp.f\
gen_G2_real.f\
global_maxi.f\
global_sumr.f\
global_sumc.f\
fmin.f\
heapsort.f\
input.f\
init_ug.f\
gaussj.f\
UxcCA.f\
UxcCA2.f\
w_line.f\
getewald.f\
getwmask.f\
getwmaskX.f\
getVrho.f\
mch_pulay.f\
mch_kerk.f\
Thomas3.f\
getNLsign.f\
getwq.f\
atomMV.f\
Etotcalc.f\
rhoIO.f\
get_ALI.f\
ugIO.f\
wqIO.f\
occup.f\
gen_Gstar_ind.f\
symmop.f\
symmopf.f\
symmcheck.f\
forcLC.f\
forcNLq.f\
forcNLr.f\
ran1.f\
dens_out.f\
system_ccfft.f\
system_csfft.f\
system_scfft.f\
system_czheev.f\
system_flush.f\
system_time.f\
system_orth_comp.f

OBJL=$(SRCL:.f=.o)

OBJL0=$(SRCL0:.f=.o)

OBJL1=$(SRCL1:.f=.o)

OBJMOD=$(MODULES:.f=.o)

$(OBJL):param.escan_real $(MODULES)

$(OBJL0):param.escan_real $(MODULES)

$(OBJL1):param.escan_real $(MODULES)

.f.o:
	$(FC) $(FFLAGS) -c $<

PEtot: $(OBJMOD) $(OBJL0) $(OBJL) $(MODULES)
	$(FC)  -o PEtot $(OBJL0) $(OBJMOD) $(OBJL) $(LIB)

Pmoment: $(OBJMOD) $(OBJL1) $(OBJL) $(MODULES) moment.f
	$(FC)  -o Pmoment $(OBJL1) $(OBJMOD) $(OBJL) $(LIB)
