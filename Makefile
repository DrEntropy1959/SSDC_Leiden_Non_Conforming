.SUFFIXES: .o .c .f90 .F90

HOSTPC:=$(shell hostname)
#=================================================================================================================
#  TAB16
#=================================================================================================================
ifeq ($(HOSTPC),tab16)
$(info ${HOSTPC})
INCLUDESDIR = -I/ump/fldmd/home/mhcarpen/OpenSourceLib/Lib-Install/include\
              -I/ump/fldmd/home/mhcarpen/OpenSourceLib/Lib-Install/petsc-3.5.2/include
FCFLAGS = -Wunused -Wmaybe-uninitialized -Wsurprising -O3 $(INCLUDESDIR)
#FCFLAGS = -Wunused -Wmaybe-uninitialized -Wsurprising -fbounds-check -O3 $(INCLUDESDIR)
#FCFLAGS = -g -Wall -Wmaybe-uninitialized -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan $(INCLUDESDIR)
CFLAGS  = -Ofast  $(INCLUDESDIR)
LFLAGS  = -L/ump/fldmd/home/mhcarpen/OpenSourceLib/Lib-Install/lib\
          -L/ump/fldmd/home/mhcarpen/OpenSourceLib/Lib-Install/petsc-3.5.2/lib
CC = gcc
FC = mpif90
LIBS = -lpetsc -lHYPRE -lumfpack -lsuperlu_dist_3.3 -lscalapack -lamd -lflapack -lfblas -lcgns -lmetis -lparmetis
#=================================================================================================================
#  Leiden
#=================================================================================================================
else ifeq ($(HOSTPC),Leiden)
$(info ${HOSTPC})
INCLUDESDIR = -I/home/carpentr/OpenSourceLib/Lib-Install/include\
              -I/home/carpentr/OpenSourceLib/Lib-Install/petsc-3.5/include
FCFLAGS = -Wunused -Wmaybe-uninitialized -Wsurprising -O3 $(INCLUDESDIR)
#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O1 -ftree-vectorizer-verbose=2 $(INCLUDESDIR)
#FCFLAGS = -g -Wall -Wmaybe-uninitialized -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)
LFLAGS = -L/home/carpentr/OpenSourceLib/Lib-Install/lib\
         -L/home/carpentr/OpenSourceLib/Lib-Install/petsc-3.5/lib
CC = gcc
FC = mpif90
LIBS = -lpetsc -lHYPRE -lumfpack -lsuperlu_dist_3.3 -lscalapack -lamd -lflapack -lfblas -lcgns -lmetis -lparmetis
#=================================================================================================================
#  dfernandez-lan
#=================================================================================================================
else ifeq ($(HOSTPC),dfernandez-lan)
$(info ${HOSTPC})
INCLUDESDIR = -I/ump/fldmd/home/ddelreyf/Documents/NASA/NEW/open_source/lib/include\
              -I/ump/fldmd/home/ddelreyf/Documents/NASA/NEW/open_source/petsc
FCFLAGS = -Wunused -Wsurprising -O1 $(INCLUDESDIR)
#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O1 -ftree-vectorizer-verbose=2 $(INCLUDESDIR)
CFLAGS = -O1  $(INCLUDESDIR)
LFLAGS = -L/ump/fldmd/home/ddelreyf/Documents/NASA/NEW/open_source/lib/lib\
         -L/ump/fldmd/home/ddelreyf/Documents/NASA/NEW/open_source/petsc
#        -L/ump/fldmd/home/ddelreyf/Documents/NASA/NEW/open_source/lib/lib\
CC = gcc
FC = mpif90
LIBS = -lpetsc -lHYPRE -lsuperlu_dist -lscalapack -lmetis -lparmetis -lcgns
else ifeq ($(HOSTPC),niavisi13-Latitude-3440)
$(info ${HOSTPC})
INCLUDESDIR = -I/home/nia-visi13/Postdoc2/open_source/deps/openmpi-1.8.3/include\
              -I/home/nia-visi13/Postdoc2/open_source/deps/petsc-3.5.2/include\
              -I/home/nia-visi13/Postdoc2/SSDC_Leiden_Non_Conforming/Libraries
FCFLAGS = -Wunused -Wmaybe-uninitialized -Wsurprising -O3 $(INCLUDESDIR)
#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O1 -ftree-vectorizer-verbose=2 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)
LFLAGS = -L/home/nia-visi13/Postdoc2/open_source/deps/openmpi-1.8.3/lib\
         -L/home/nia-visi13/Postdoc2/open_source/deps/petsc-3.5.2/lib\
         -L/home/nia-visi13/Postdoc2/SSDC_Leiden_Non_Conforming/Libraries
CC = gcc
FC = mpif90
LIBS = -lpetsc -lHYPRE -lumfpack -lsuperlu_dist_3.3 -lscalapack -lamd -lflapack -lfblas -lcgns -lmetis -lparmetis
#=================================================================================================================
#  computer-dcdelrey
#=================================================================================================================
else ifeq ($(HOSTPC),computer-dcdelrey)
$(info ${HOSTPC})
INCLUDESDIR = -I/home/dcdelrey/Documents/Postdoc2/open_source/deps/openmpi-1.8.3/include\
              -I/home/dcdelrey/Documents/Postdoc2/open_source/deps/petsc-3.5.2/include\
              -I/home/dcdelrey/Documents/Postdoc2/SSDC_Leiden_Non_Conforming/Libraries
FCFLAGS = -Wunused -Wmaybe-uninitialized -Wsurprising -O3 $(INCLUDESDIR)
#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O1 -ftree-vectorizer-verbose=2 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)
LFLAGS = -L/home/dcdelrey/Documents/Postdoc2/open_source/deps/openmpi-1.8.3/lib\
         -L/home/dcdelrey/Documents/Postdoc2/open_source/deps/petsc-3.5.2/lib\
         -L/home/dcdelrey/Documents/Postdoc2/SSDC_Leiden_Non_Conforming/Libraries
CC = gcc
FC = mpif90
LIBS = -lpetsc -lHYPRE -lumfpack -lsuperlu_dist_3.3 -lscalapack -lamd -lflapack -lfblas -lcgns -lmetis -lparmetis
else
$(info Hostname Not Found)
endif

SRCS = precision_vars.f90\
       datatypes.f90\
       variables.f90\
       referencevariables.f90\
       controlvariables.f90\
       nsereferencevariables.f90\
       collocationvariables.f90\
       SSWENOvariables.f90\
       petscvariables.F90\
       initcollocation.f90\
       fileio.f90\
       mpimod.F90\
       CSRlocalvariables.f90\
       unary_mod.f90\
       tools.f90\
       interpolation.f90\
       jacobian_matrix_implicit_ts_variables.f90\
       matvec_module.f90\
       initialize_CSR.f90\
       tools_IO.f90\
       time_average.f90\
       SSWENO_routines.f90\
       restart_simulation.f90\
       eispack_module.f90\
       initgrid.F90\
       navierstokes.f90\
       time_integ_coeff.f90\
       physics_driver.f90\
       errorestimation.f90\
       dkinetic_energy_dt_enstrophy.f90\
       jacobian_matrix_implicit_ts.f90\
       implicit_residual.f90\
       petsc_snes_solver.F90\
       binary_sizes.f90\
       write_solution_file.f90\
       aerodynamic_coefficients.f90\
       error_bc_no_slip_wall.f90\
       error_heat_entropy_flow_wall_bc.f90\
       timeinteg.f90\
       non_conforming.f90\
       polyinit.f90\
       physicsindependent.f90\
       ssdcsolver.f90

#OBJS1 = $(SRCS:.f90=.o)
#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = precision_vars.o\
       datatypes.o\
       variables.o\
       referencevariables.o\
       controlvariables.o\
       nsereferencevariables.o\
       collocationvariables.o\
       SSWENOvariables.o\
       petscvariables.o\
       initcollocation.o\
       fileio.o\
       mpimod.o\
       CSRlocalvariables.o\
       unary_mod.o\
       tools.o\
       interpolation.o\
       jacobian_matrix_implicit_ts_variables.o\
       matvec_module.o\
       initialize_CSR.o\
       tools_IO.o\
       time_average.o\
       SSWENO_routines.o\
       restart_simulation.o\
       eispack_module.o\
       initgrid.o\
       navierstokes.o\
       time_integ_coeff.o\
       physics_driver.o\
       errorestimation.o\
       dkinetic_energy_dt_enstrophy.o\
       jacobian_matrix_implicit_ts.o\
       implicit_residual.o\
       petsc_snes_solver.o\
       binary_sizes.o\
       write_solution_file.o\
       aerodynamic_coefficients.o\
       error_bc_no_slip_wall.o\
       error_heat_entropy_flow_wall_bc.o\
       timeinteg.o\
       non_conforming.o\
       polyinit.o\
       physicsindependent.o\
       ssdcsolver.o\
       metiselementpartition.o

#LINK = gfortran

TARGET = SSDCNSE

.PHONY: depend clean

#metiselementpartition.o: metiselementpartition.c


all: $(TARGET)
	@echo  SSDC code has been compiled, creating executable named SSDCNSE


$(TARGET) : $(OBJS)
	$(FC) $(FCFLAGS) $(INCLUDESDIR) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) -c $< -o $@

%.o : %.F90
	$(FC) -c $(FCFLAGS) -c $< -o $@

%.o : %.c
	$(CC) -c $(CFLAGS) -c $< -o $@


clean:
	rm -f *.o *.mod *~ $(TARGET)

depend: $(SRCS)
	makedepend $(INCLUDESDIR) $^

# DO NOT DELETE THIS LINE -- make depend needs it				
