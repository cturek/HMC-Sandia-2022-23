#
##
## Makefile for the serial SeqQuest code
##
## Last revised: 11Feb12-PAS: for version 2.62 (clean-up)
## Configured for serial: g77
##
#

######
###### Choose compiler (pick one and only one):
######
## general:
# FC = f77
## g77
#FC = g77
## pgi compiler:
#FC = /opt/pgi/linux86-64/7.2-1/bin/pgf77
## mpi:
# FC = mpif77
FC = gfortran
 
######
###### Set compiler flags (pick one and only one):
######
## general/preferred:
# FFLAGS = -O -fast
## general/when -fast is not available, or with pgi compiler:
# FFLAGS = O
FFLAGS = -O2 -std=legacy
## for some g77 type-change (type-casting?) causes compile failure:
# FFLAGS = -O -fno-globals
## on SGI machines:
# FFLAGS = -O3 -mips4 -OPT:IEEE_arithmetic=3:alias=restrict


LFLAGS =
LIBS = 

######
###### Math libraries (pick one and only one):
######
#
# Select appropriate math library to link in.
# Internally, code default is to make use of blas/3 routines,
# an option that is flag-controlled (doblas3).  With an
# optimized blas library, this gives the best performance,
# but without an optimized library, using blas/3 routines
# gives much slower execution.  If your platform does not
# have an optimized library and you must use the explicit
# source provided below, edit subroutine "dat0opt" in the
# file lcaosubs.f to set default "doblas3" to .false. for
# most efficient execution.
#
## explicit source as provided:
 MATH = lapack.o blas.o

#########################################################################
#########################################################################
#########################################################################

######
###### Source files
######
## main program
MVRS=lcaomain.o
#MVRSTP=lcaomain_tp.o
## program subroutines
SUBS=lcaosubs.o
#TPSUBS=lcaosubs_tp.o
#SUBS=subs/subs.a
#TPSUBS=subs_tp/subs_tp.a

######
###### Platform-dependent configuration file (serial or mpi-parallel)
######
## general serial code
SUBSU = utl_gen.o
## general MPI code
#SUBSU = utl_mpi.o

######
###### Where it all actually happens
######

#seqquest-tp: quest_ompi.x

seqquest: lcao.x


clean:
	rm -f *.o lcao.x lcao_m.x

######
###### The actual builds
######

lcao.x : ${MVRS} ${SUBS} ${SUBSU} ${MATH}
	$(FC) $(LFLAGS) -o $@  ${MVRS} ${SUBS} ${SUBSU} ${MATH}

lcao_m.x : ${MVRS} ${SUBS} ${SUBSU} lapack.o blas.o
	$(FC) $(LFLAGS) -o $@  $^
