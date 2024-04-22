#########################################################################
SHELL = /bin/bash
FFLAG = -O5 -w
# FFTIDIR  = -I/THL6/home/hegw/Shibeiji/FLASH/FlashLib/FFTW3/Double/include
# FFTLDIR  = -L/THL6/home/hegw/Shibeiji/FLASH/FlashLib/FFTW3/Double/lib
IDIR  = -I/home/wuct/lapack-3.12.0/LAPACKE/include
LDIR  = -llapack -lrefblas -lz -lm
FCOMP =  mpif90 -c -g ${FFLAG} ${IDIR} ${FFTIDIR}
LINK  =  mpif90
LIBS  = -lfftw3 -lm

#module name
MOD   = global

OBJ   =  main.o

.SUFFIXES: .o .f90 .mod

.f90.o:
	${FCOMP} $*.f90

main:   ${OBJ}
	${LINK} ${OBJ} ${LDIR} ${FFTLDIR} ${LIBS} -o main


cleanall:
	rm -f *.o
	rm -rf *.mod
	rm -f *-pp.f90
	rm -rf main