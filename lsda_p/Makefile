#Makefile for PEtot on IBM SP

FC=mpxlf90_r
LDR=mpxlf90_r
#FFLAGS = -O2 -qfixed -qarch=pwr3 -I$$SCRATCH -qmoddir=$$SCRATCH
FFLAGS = -O2  -qfixed -qarch=pwr3  -qmaxmem=4096
#FFLAGS = -g -pg -qfixed -qarch=pwr3 -qmaxmem=4096
#FFLAGS = -g -qfixed -qarch=pwr3 -I$$SCRATCH -qmoddir=$$SCRATCH
LIB = $$LAPACK -lessl


MODULES=load_data.f \
fft_data.f \
data.f

SRCL0=mainMV.f

SRCL1=moment.f

SRCL=Hpsi_comp.f\
CG_comp.f\
CG_new.f\
getpot2L.f\
getpot3L.f\
getpot4L.f\
getpot4_force.f\
getpot5L.f\
getpot5_force.f\
GGAPBE.f\
d3fft_comp.f\
cfft.f\
cfftd.f\
diag_comp.f\
djacobi.f\
convert_SLvr.f\
d3fft_real2.f\
d3fft_real2L.f\
d3fft_real2L2.f\
fftprep_comp.f\
fftprep_real2.f\
fftprep_real2L.f\
fftprep_real2L2.f\
fwdcpfft_comp.f\
fwdcpfft2.f\
fwdcpfft2L.f\
fwdcpfft2L2.f\
invcpfft_comp.f\
invcpfft2.f\
invcpfft2L.f\
invcpfft2L2.f\
gen_G_comp.f\
gen_G2_real.f\
gen_G2L_real.f\
gen_G2L2_real.f\
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
getwmask_q.f\
add_rho_beta.f\
getVrhoL.f\
mch_pulay.f\
mch_kerk.f\
Thomas3.f\
getNLsign.f\
getwq.f\
atomMV.f\
Etotcalc.f\
rhoIO.f\
ugIO.f\
wqIO.f\
beta_psiIO.f\
get_ALI.f\
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
densWr_out.f\
readusp_head.f\
w_line_vwr.f\
w_line_usp.f\
LegendreSP.f\
clebsch_gordan.f\
get_Dij.f\
system_orth_comp.f\
system_ccfft.f\
system_csfft.f\
system_scfft.f\
system_czheev.f\
system_flush.f\
getwmask_dq.f\
get_VdqdR.f\
getvcoul.f\
convert_2LtoL.f\
getV_Hartree.f\
getewald3D.f\
getewald2D.f\
forcLC2.f\
getEextV.f\
getrho_only.f\
write_wg.f\
system_time.f

OBJL=$(SRCL:.f=.o)

OBJL0=$(SRCL0:.f=.o)

OBJL1=$(SRCL1:.f=.o)

OBJMOD=$(MODULES:.f=.o)

$(OBJL):param.escan_real $(MODULES)

$(OBJL0):param.escan_real $(MODULES)

$(OBJL1):param.escan_real $(MODULES)

.f.o:
	$(FC) $(FFLAGS) -c $<

PEtot_U: $(OBJMOD) $(OBJL0) $(OBJL) $(MODULES)
	$(LDR) -bmaxdata:0x80000000  -o PEtot_U $(OBJL0) $(OBJMOD) $(OBJL) $(LIB)

Pmoment: $(OBJMOD) $(OBJL1) $(OBJL) $(MODULES) moment.f
	$(LDR) -bmaxdata:0x80000000 -o Pmoment $(OBJL1) $(OBJMOD) $(OBJL) $(LIB)

clean:
	rm -f *.o
