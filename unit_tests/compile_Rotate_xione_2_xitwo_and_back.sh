#!/bin/sh
# compile Rotate_xione_2_xitwo_and_back.f90  
# go up one directory
cd ..

SRC=$PWD # where are the source files
OBJ=$PWD # where are the mod files
# -O0 = no optimization, -Og = optimization for debugging
#-Wall = enable all warnings
#-fimplicit-none = no implicit typing allowed
#-fcheck=all = enable run-time tests, such as array bounds check
#-fbacktrace = if a deadly signal is encountered output a backtrace of the error
FLAGS="-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan"
#BIN      # where to put the program (not currently used)

#-- go back to unit_tests directory
cd unit_tests

#compile the mod files
gfortran -J$OBJ $FLAGS -c $SRC/precision_vars.f90 -o $OBJ/precision_vars.o
gfortran -J$OBJ $FLAGS -c $SRC/non_conforming.f90 -o $OBJ/non_conforming.o
gfortran -J$OBJ $FLAGS -c $SRC/referencevariables.f90 -o $OBJ/referencevariables.o
gfortran -J$OBJ $FLAGS -c $SRC/initcollocation.f90 -o $OBJ/initcollocation.o
gfortran -J$OBJ $FLAGS -c $SRC/collocationvariables.f90 -o $OBJ/collocationvariables.o

#compile the program
gfortran $FLAGS -I$OBJ -c Rotate_xione_2_xitwo_and_back.f90
# make the executable
gfortran -o Rotate_xione_2_xitwo_and_back \
Rotate_xione_2_xitwo_and_back.o \
$OBJ/precision_vars.o \
$OBJ/non_conforming.o \
$OBJ/referencevariables.o \
$OBJ/initcollocation.o \
$OBJ/collocationvariables.o

