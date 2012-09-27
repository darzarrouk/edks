FCOMP  = gfortran44
#FCOMP = ifort
BINDIR = ../../bin
#FFLAGS = -g90 -C
#FFLAGS = -O4 -C
#FFLAGS = -g -bytereclen
#FFLAGS = -Ofast -bytereclen
#FFLAGS = -O3 -bytereclen
#FFLAGS = -04
# uncomment for use with gfortran
FFLAGS = -frecord-marker=4 -O4
# uncomment for use with ifort
#FFLAGS = ''


all:		tab4E build_edks

tab4Eg:           tab4E.o ken.o ertudpsv.o calcqwvE.o calc_coeffsE.o \
	   	nrtype.o nr.o nrutil.o bessj0.o bessj1.o bessj.o \
		calc_BCD0E.o discotabE.o read_par_mod.o acotabE.o
		${FCOMP}  $(FFLAGS) -o  ${BINDIR}/tab4Eg tab4E.o ken.o ertudpsv.o calcqwvE.o calc_coeffsE.o nrtype.o nr.o nrutil.o bessj0.o bessj1.o bessj.o calc_BCD0E.o discotabE.o read_par_mod.o acotabE.o

tab4E:           modules.o tab4E.o ken.o ertudpsv.o calcqwvE.o calc_coeffsE.o \
	   	nrtype.o nr.o nrutil.o bessj0.o bessj1.o bessj.o \
		calc_BCD0E.o discotabE.o read_par_mod.o acotabE.o
		${FCOMP}  $(FFLAGS) -o  ${BINDIR}/tab4E modules.o tab4E.o ken.o ertudpsv.o calcqwvE.o calc_coeffsE.o \
				   nrtype.o nr.o nrutil.o bessj0.o bessj1.o bessj.o \
				   calc_BCD0E.o discotabE.o read_par_mod.o acotabE.o

modules.o:		               modules.f90
		${FCOMP}  $(FFLAGS) -c modules.f90

tab4E.o:		               tab4E.f90
		${FCOMP}  $(FFLAGS) -c tab4E.f90

discotabE.o:	                       discotabE.f90
		${FCOMP}  $(FFLAGS) -c discotabE.f90

read_par_mod.o:	                       read_par_mod.f90
		${FCOMP}  $(FFLAGS) -c read_par_mod.f90

acotabE.o:	                       acotabE.f90
		${FCOMP}  $(FFLAGS) -c acotabE.f90

ken.o:		                       ken.f90
		${FCOMP}  $(FFLAGS) -c ken.f90

ertudpsv.o:	                       ertudpsv.f90
		${FCOMP}  $(FFLAGS) -c ertudpsv.f90

calcqwvE.o:	                       calcqwvE.f90
		${FCOMP}  $(FFLAGS) -c calcqwvE.f90

mecan.o:	                       mecan.f90
		${FCOMP}  $(FFLAGS) -c mecan.f90

mecamom.o:	                       mecamom.f90
		${FCOMP}  $(FFLAGS) -c mecamom.f90

calc_BCD0E.o:	                       calc_BCD0E.f90
		${FCOMP}  $(FFLAGS) -c calc_BCD0E.f90

calc_coeffsE.o:	                       nrtype.f90 nr.f90 nrutil.f90 calc_coeffsE.f90
		${FCOMP}  $(FFLAGS) -c nrtype.f90 nr.f90 nrutil.f90 calc_coeffsE.f90

nrutil.o:	                       nrutil.f90
		${FCOMP}  $(FFLAGS) -c nrutil.f90

nrtype.o:	                       nrtype.f90
		${FCOMP}  $(FFLAGS) -c nrtype.f90

nr.o:		                       nr.f90
		${FCOMP}  $(FFLAGS) -c nr.f90

bessj0.o:                   nrtype.f90 nr.f90 nrutil.f90 bessj0.f90 
		${FCOMP}  $(FFLAGS) -c nrtype.f90 nr.f90 nrutil.f90 bessj0.f90
bessj1.o:                   nrtype.f90 nr.f90 nrutil.f90 bessj1.f90 
		${FCOMP}  $(FFLAGS) -c nrtype.f90 nr.f90 nrutil.f90 bessj1.f90

bessj.o:                    nrtype.f90 nr.f90 nrutil.f90 bessj.f90 
		${FCOMP}  $(FFLAGS) -c nrtype.f90 nr.f90 nrutil.f90 bessj.f90

clean:
		\rm *.o *.mod

build_edks:		                                    build_edks.o write_edks.o read_par_mod.o
		${FCOMP}  $(FFLAGS) -o ../../bin/build_edks build_edks.o write_edks.o read_par_mod.o

build_edks.o:				build_edks.f90
		 ${FCOMP}  $(FFLAGS) -c build_edks.f90

write_edks.o:				write_edks.f90
		 ${FCOMP}  $(FFLAGS) -c write_edks.f90
