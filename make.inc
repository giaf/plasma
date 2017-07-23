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

#CFLAGS   = ${FPIC} -O2 -std=c99 -fopenmp -DHAVE_OPENMP_DEPEND -DHAVE_OPENBLAS -I/opt/openblas/include -I/opt/blasfeo/include -DBLASFEO
CFLAGS   = ${FPIC} -O2 -std=c99 -fopenmp -DHAVE_OPENMP_DEPEND -DHAVE_MKL -m64 -I/opt/intel/mkl/include -I/opt/blasfeo/include #-DBLASFEO
FCFLAGS  = ${FPIC} -O2 -std=f2008 -fopenmp
LDFLAGS  = ${FPIC}  -fopenmp
#LIBS     = /opt/openblas/lib/libopenblas.a /opt/blasfeo/lib/libblasfeo.a -lm
LIBS     = -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a -Wl,--end-group -ldl -lpthread /opt/blasfeo/lib/libblasfeo.a -lm

# Enable or disable compiling Fortran 2008 interfaces into PLASMA library
# 0 - disabled
# 1 - enabled; build Fortran interfaces and examples
fortran ?= 0

# where to install PLASMA
prefix ?= /opt/plasma

# one of: aix bsd c89 freebsd generic linux macosx mingw posix solaris
# usually generic is fine
lua_platform ?= generic
