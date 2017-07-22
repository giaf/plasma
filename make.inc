# PLASMA make.inc template, processed by configure.py
#
# PLASMA is a software package provided by:
# University of Tennessee, US,
# University of Manchester, UK.

# ------------------------------------------------------------------------------
# programs and flags

CC       = gcc-6
FC       = gfortran
RANLIB   = ranlib
AR       = ar

# Use -fPIC to make shared (.so) and static (.a) libraries;
# can be commented out if making only static libraries.
FPIC     = -fPIC

CFLAGS   = ${FPIC} -std=c99 -fopenmp -DHAVE_OPENMP_DEPEND -DHAVE_OPENMP_DEPEND -DHAVE_OPENBLAS -I/opt/openblas/include
FCFLAGS  = ${FPIC} -std=f2008 -fopenmp
LDFLAGS  = ${FPIC}  -fopenmp
LIBS     = /opt/openblas/lib/libopenblas.a -lm

# Enable or disable compiling Fortran 2008 interfaces into PLASMA library
# 0 - disabled
# 1 - enabled; build Fortran interfaces and examples
fortran ?= 0

# where to install PLASMA
prefix ?= /opt/plasma

# one of: aix bsd c89 freebsd generic linux macosx mingw posix solaris
# usually generic is fine
lua_platform ?= generic
