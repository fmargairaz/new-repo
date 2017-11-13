######################################################################################
#
# Author: Fabien Margairaz
#
######################################################################################

#define shell
SHELL = /bin/sh

#architecture ---------------------------- please modify
ARCH = KINGSPEAK_INTEL

#code options
OPTIONS = #-DDNS_MOD

#dealiasing scheme
DSLG = PADD
#DSLG = TRUNC
#DSLG = FSMTHG
#DSLG=NODLSG

ifeq ($(DSLG),PADD)
  #executable name
  EXE = les_code_mpi.x
  OPTIONS += -DDEA_PADD
endif
ifeq ($(DSLG),TRUNC)
  EXE = les_code_mpi_trunc.x
  OPTIONS += -DDEA_TRUNC
endif
ifeq ($(DSLG),NODLSG)
  EXE = les_code_mpi_nodlsg.x
  OPTIONS += -DDEA_NODLSG
endif
ifeq ($(DSLG),FSMTHG)
  EXE = les_code_mpi_Fsmthg.x
  OPTIONS += -DDEA_FSMTHG
endif


#directory for the .o files
OBJDIR = ./.obj
MODDIR = ./.mod
EXEDIR = ..


######################################################################################
#  setting compilers and libraries depending on the pc
######################################################################################

#settings for BELLATRIX with INTEL compiler (mpiifort)
ifeq ($(ARCH),BELLATRIX_INTEL)
  FC = mpiifort
  FFLAGS = -fpp $(FPPF) -mcmodel medium -openmp -module $(MODDIR) -I$(OBJDIR) #-openmp
  LDFLAGS = -mcmodel=medium

  FFTPATH = /home/giometto/libraries
  DECPATH = /home/giometto/fft_pjt_fmargairaz/2decomp_fft
  VTKPATH = /home/giometto/fft_pjt_fmargairaz/lib_vtk
  FUTPATH = /home/giometto/fft_pjt_fmargairaz/futils
  HDFPATH = /home/giometto/fft_pjt_fmargairaz/hdf5

  LIBS = -L$(FFTPATH) -lfftw3 -L$(DECPATH)/lib -l2decomp_fft -L$(VTKPATH) -llib_vtk_io 
  LIBS += -L$(FUTPATH)/lib -lfutils -L$(HDFPATH)/lib -lhdf5_fortran -lhdf5 -lz
  INCS = -I$(FFTPATH)/include -I$(DECPATH)/include -I$(VTKPATH)/ 
  INCS += -I$(FUTPATH)/include -I$(HDFPATH)/include

  #dependencies
  MAKEDEP = /home/giometto/fft_pjt_fmargairaz/makedepf90
endif

#settings for KINGSPEAK with INTEL compiler (mpiifort)
ifeq ($(ARCH),KINGSPEAK_INTEL)
  FC = mpiifort
  #FC = mpif90	 
  #FFLAGS = -cpp -module $(MODDIR) -I$(OBJDIR) -traceback -check bound -mcmodel=medium -lhugetlbfs
  FFLAGS = -cpp -O3 -mcmodel=medium -lhugetlbfs -module $(MODDIR) -I$(OBJDIR) -qopenmp 
  #LDFLAGS = -mcmodel=medium

  LIBPATH = /uufs/chpc.utah.edu/common/home/u0917142/libraries
  FFTPATH = /uufs/chpc.utah.edu/sys/pkg/fftw/std_intel

  LIBS = -L$(LIBPATH)/lib -l2decomp_fft -L$(FFTPATH)/lib -lfftw3
  INCS = -I$(LIBPATH)/include -I$(FFTPATH)/include 

  #dependencies
  MAKEDEP = /uufs/chpc.utah.edu/common/home/u0917142/bin/makedepf90

  #LIBS += -L$(FUTPATH)/lib -lfutils -L$(HDFPATH)/lib -lhdf5_fortran -lhdf5 -lz
  #INCS += -I$(FUTPATH)/include -I$(HDFPATH)/include


endif


#################################################################################
#  various checks and launch the makefile for compilation
#################################################################################

# preprocessor flags settings
COMPSTR = '$$(FC) $$(OPTIONS) -c $$(INCS) -o $$@ $$(FFLAGS) $$<'

#include the dependency list created by makedepf90 below
include .depend

.depend:
	mkdir $(OBJDIR) 
	mkdir $(MODDIR)	
	$(MAKEDEP) -r $(COMPSTR) -b $(OBJDIR) -o $(EXE) *.f90 > .depend


#################################################################################
# cleaning and extra useful commands
#################################################################################

.PHONY : clean

#cleaning out everything
clean:
	rm .depend
	rm -r $(OBJDIR)
	rm -r $(MODDIR)
	rm -v *.x
	@echo "All compilation files wiped out!"
