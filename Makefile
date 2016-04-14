#
# --- HYCOM 2.2 makefile 
#
# --- Stand-alone HYCOM, or HYCOM ESMF component, or HYCOM+CICE.
#
# --- Tunable parameters in ../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
.f:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

#include ../config/$(ARCH)_$(TYPE)
include Aintelsse-impi-sm-nc-btrmas_mpi

MODS =   mod_dimensions.o mod_xc.o mod_za.o mod_cb_arrays.o mod_pipe.o \
         mod_incupd.o \
         mod_floats.o mod_stokes.o mod_tides.o mod_nc.o mod_mean.o mod_archiv.o \
         mod_tsadvc.o mod_momtum.o \
         mod_hycom.o

MODD =   mod_dimensions.o mod_xc.o mod_za.o mod_cb_arrays.o mod_pipe.o \
         mod_incupd.o \
         mod_floats.o mod_stokes.o mod_tides.o mod_nc.o mod_mean.o mod_archiv.o \
         mod_tsadvc.o mod_momtum.o \
         mod_hycom_dummy.o

OBJS =	         barotp.o  bigrid.o blkdat.o  cnuity.o convec.o \
	diapfl.o dpthuv.o  dpudpv.o forfun.o  geopar.o hybgen.o \
	icloan.o inicon.o inigiss.o inikpp.o   inimy.o latbdy.o \
	matinv.o mxkprf.o  mxkrt.o  mxkrtm.o  mxpwp.o \
	overtn.o poflat.o  prtmsk.o  psmoo.o restart.o \
	thermf.o trcupd.o  \
       machine.o  wtime.o machi_c.o  isnan.o sgefs.o asselin.o latbdt_river.o

hycom:	$(MODS) $(OBJS) hycom.o
	$(LD)  $(LDFLAGS) -o hycom  hycom.o $(MODS) $(OBJS) $(EXTRALIBS)

esmf:	$(MODS) $(OBJS)
	@echo "--- ESMF hycom component has been built ---"

hycom_cice:	$(MODS) $(OBJS) mod_OICPL.o hycom_cice.o
	$(LD)  $(LDFLAGS) -o hycom_cice \
                             hycom_cice.o mod_OICPL.o \
                             $(MODS) $(OBJS) \
                             ${CICE_DIR}/esmf/compile/*.o \
                             $(EXTRALIBS)

dummy_cice:	$(MODS) $(OBJS) mod_OICPL.o dummy_cice.o
	$(LD)  $(LDFLAGS) -o dummy_cice \
                             dummy_cice.o mod_OICPL.o \
                             $(MODD) $(OBJS) \
                             ${CICE_DIR}/esmf/compile/*.o \
                             $(EXTRALIBS)

hycom.o:        hycom.F       mod_hycom.o
hycom_cice.o:   hycom_cice.F  mod_hycom.o       mod_OICPL.o
dummy_cice.o:   dummy_cice.F  mod_hycom_dummy.o mod_OICPL.o

asselin.o:   asselin.F        mod_cb_arrays.o stmt_fns.h
barotp.o:  barotp.F  mod_xc.o mod_cb_arrays.o                     mod_pipe.o \
	                                                          mod_tides.o \
	                                                          mod_stokes.o
bigrid.o:  bigrid.f  mod_xc.o 
blkdat.o:  blkdat.F  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_incupd.o \
	                                                          mod_floats.o \
	                                                          mod_tides.o \
	                                                          mod_stokes.o
cnuity.o:  cnuity.F  mod_xc.o mod_cb_arrays.o                     mod_pipe.o \
	                                                          mod_stokes.o
convec.o:  convec.F  mod_xc.o mod_cb_arrays.o stmt_fns.h
diapfl.o:  diapfl.F  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_stokes.o
dpthuv.o:  dpthuv.F  mod_xc.o mod_cb_arrays.o
dpudpv.o:  dpudpv.F  mod_xc.o 
forfun.o:  forfun.f  mod_xc.o mod_cb_arrays.o            mod_za.o
geopar.o:  geopar.F  mod_xc.o mod_cb_arrays.o stmt_fns.h mod_za.o
hybgen.o:  hybgen.F  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o
icloan.o:  icloan.F  mod_xc.o mod_cb_arrays.o stmt_fns.h
inicon.o:  inicon.F  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o
inigiss.o: inigiss.f mod_xc.o mod_cb_arrays.o stmt_fns.h
inikpp.o:  inikpp.f  mod_xc.o mod_cb_arrays.o stmt_fns.h
inimy.o:   inimy.f   mod_xc.o mod_cb_arrays.o stmt_fns.h
isnan.o:   isnan.F
latbdy.o:  latbdy.F  mod_xc.o mod_cb_arrays.o                     mod_tides.o
latbdt_river.o:  latbdt_river.F  mod_xc.o mod_cb_arrays.o         mod_nc.o
machine.o: machine.F
machi_c.o: machi_c.c
matinv.o:  matinv.f  mod_xc.o mod_cb_arrays.o
mxkprf.o:  mxkprf.F  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
	                                                          mod_stokes.o
mxkrt.o:   mxkrt.F   mod_xc.o mod_cb_arrays.o stmt_fns.h
mxkrtm.o:  mxkrtm.F  mod_xc.o mod_cb_arrays.o stmt_fns.h
mxpwp.o:   mxpwp.F   mod_xc.o mod_cb_arrays.o stmt_fns.h
overtn.o:  overtn.F  mod_xc.o mod_cb_arrays.o
poflat.o:  poflat.f
prtmsk.o:  prtmsk.f
psmoo.o:   psmoo.f   mod_xc.o 
restart.o: restart.f mod_xc.o mod_cb_arrays.o            mod_za.o mod_tides.o
sgefs.o:   sgefs.f
thermf.o:  thermf.F  mod_xc.o mod_cb_arrays.o stmt_fns.h
trcupd.o:  trcupd.F  mod_xc.o mod_cb_arrays.o                     mod_pipe.o
wtime.o:   wtime.F
mod_hycom.o: \
        mod_hycom.F  mod_xc.o mod_cb_arrays.o            mod_za.o mod_pipe.o \
	                                                          mod_incupd.o \
	                                                          mod_mean.o \
	                                                          mod_floats.o \
	                                                          mod_momtum.o \
	                                                          mod_tsadvc.o \
	                                                          mod_stokes.o
mod_hycom_dummy.o: \
        mod_hycom_dummy.F  mod_xc.o mod_cb_arrays.o      mod_za.o mod_pipe.o \
	                                                          mod_incupd.o \
	                                                          mod_mean.o \
	                                                          mod_floats.o \
	                                                          mod_momtum.o \
	                                                          mod_tsadvc.o
mod_cb_arrays.o: \
        mod_cb_arrays.F mod_dimensions.o
mod_momtum.o: \
	mod_momtum.F mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
	                                                          mod_tides.o \
	                                                          mod_stokes.o
mod_tsadvc.o: \
	mod_tsadvc.F mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o
mod_incupd.o: \
        mod_incupd.F mod_xc.o mod_cb_arrays.o            mod_za.o
mod_floats.o: \
        mod_floats.F mod_xc.o mod_cb_arrays.o            mod_za.o mod_pipe.o \
	                                                          mod_stokes.o
mod_pipe.o: \
        mod_pipe.F   mod_xc.o mod_cb_arrays.o                     mod_stokes.o
mod_stokes.o: \
        mod_stokes.F mod_xc.o mod_cb_arrays.o            mod_za.o
mod_tides.o: \
        mod_tides.F  mod_xc.o mod_cb_arrays.o            mod_za.o
mod_mean.o: \
        mod_mean.F   mod_xc.o mod_cb_arrays.o            mod_za.o mod_stokes.o
mod_archiv.o: \
        mod_archiv.F mod_xc.o mod_cb_arrays.o            mod_za.o mod_stokes.o
mod_nc.o: \
        mod_nc.F     mod_xc.o mod_cb_arrays.o            mod_za.o 

mod_dimensions.o:   mod_dimensions.F dimensions.h dimensions_relo.h
mod_xc.o: mod_xc.F  mod_dimensions.o mod_xc_sm.h mod_xc_mp.h
mod_za.o: mod_za.F  mod_xc.o         mod_za_sm.h mod_za_mp.h mod_za_mp1.h mod_za_zt.h

mod_OICPL.o: mod_OICPL.F
