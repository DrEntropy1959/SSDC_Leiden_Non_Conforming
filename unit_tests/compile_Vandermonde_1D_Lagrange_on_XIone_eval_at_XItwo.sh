#!/bin/sh
# compile Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo.f90 
# go up one directory
cd ..

SRC=$PWD # where are the source files
OBJ=$PWD # where are the mod files
#BIN      # where to put the program (not currently used)

#-- go back to unit_tests directory
cd unit_tests

#compile the mod files
gfortran -J$OBJ -c $SRC/precision_vars.f90 -o $OBJ/precision_vars.o
gfortran -J$OBJ -c $SRC/non_conforming.f90 -o $OBJ/non_conforming.o
#compile the program
gfortran -I$OBJ -c Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo.f90 -o \ 
Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo.o
# make the executable
gfortran -o Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo \ 
Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo.o \
$OBJ/precision_vars.o \
$OBJ/non_conforming.o

