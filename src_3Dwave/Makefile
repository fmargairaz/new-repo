######################################################################################
#
# MAKEFILE WITH AUTOMATIC DEPENDENCIES for program prova.f90
#
######################################################################################
#
# Date of last update: 07/01/2013
#
#
# Author: Marco Giometto
#
# Environment: Ubuntu 12.04
#
######################################################################################

#define shell
SHELL = /bin/sh

#architecture ---------------------------- please modify
#ARCH = MAC_GNU
ARCH = KINGSPEAK_INTEL

#executable name
EXE = test_code_mpi.x

#directory for the .o files
OBJDIR = ./.obj
MODDIR = ./.mod
EXEDIR = ..

######################################################################################
#  setting compilers and libraries depending on the pc
######################################################################################


#settings for MAC with GNU compiler (gfortran)
ifeq ($(ARCH),MAC_GNU)
  #FC = gfortran
  FC = mpif90
  FFLAGS = -O3 -cpp -m64 -Wall -fbacktrace -g -fbounds-check -mcmodel=medium -I$(OBJDIR) -J$(MODDIR) 
  FFTPATH = /usr/local
  INCPATH = /usr/local
  DECPATH = /users/fabienmargairaz/software/2decomp_fft
  LIBS = -L$(fftpath) -lfftw3 -lfftw3f -L$(decpath)/lib -l2decomp_fft 
  INC = -I$(fftpath)/include -I$(decpath)/include
  #dependencies
  MAKEDEP = makedepf90
endif

#settings for PC with INTEL compiler (ifort)
ifeq ($(ARCH),PC_INTEL)
  FC = ifort
  FFLAGS = -O3 -traceback -fpp $(FPPF) -mcmodel medium -module $(MODDIR) -I$(OBJDIR) #-openmp
  LDFLAGS = -mcmodel medium
  LIBPATH = -L/home/mg/Phd/libraries/bin_gfortran
  LIBS = $(LIBPATH) -lfftw3 
endif

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

#settings for KINGSPEAK with INTEL compiler (mpif90)
ifeq ($(ARCH),KINGSPEAK_INTEL)
  FC = mpif90
  FFLAGS = -mcmodel medium -openmp -module $(MODDIR) -I$(OBJDIR) #-openmp
  LDFLAGS = -mcmodel=medium

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
COMPSTR = '$$(FC) -c $$(INCS) -o $$@ $$(FFLAGS) $$<'

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
	rm $(EXE)
	@echo "All compilation files wiped out!"
