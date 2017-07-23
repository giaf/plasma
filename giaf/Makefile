CC = gcc
#CFLAGS = -O2 -I../test/ -I/opt/plasma/include -I/opt/lua/include -I/opt/openblas/include -fopenmp
CFLAGS = -O2 -I../test/ -I/opt/plasma/include -I/opt/lua/include -m64 -I/opt/intel/mkl/include -fopenmp

#LIBS = ../lib/libplasma.a ../lib/libcoreblas.a /opt/lua/lib/liblua.a /opt/openblas/lib/libopenblas.a /opt/blasfeo/lib/libblasfeo.a -lm -fopenmp
LIBS = ../lib/libplasma.a ../lib/libcoreblas.a /opt/lua/lib/liblua.a -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a -Wl,--end-group -ldl -lpthread /opt/blasfeo/lib/libblasfeo.a -lm -fopenmp
#LIBS = ../lib/libplasma.a ../lib/libcoreblas.a /opt/lua/lib/liblua.a -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread /opt/blasfeo/lib/libblasfeo.a -lm -fopenmp

#OBJS = test_dgemm.o
#OBJS = test_dpotrf.o
OBJS = test_dpotrf_nopack.o

all: $(OBJS)

test:
	$(CC) -o test.out $(OBJS) $(LIBS)

run:
	./test.out

clean:
	rm -r *.o
	rm -r test.out
