CC = gcc
CFLAGS = -O2 -I../test/ -I/opt/plasma/include -I/opt/lua/include -I/opt/openblas/include

LIBS = /opt/plasma/lib/libplasma.a /opt/plasma/lib/libcoreblas.a /opt/lua/lib/liblua.a /opt/openblas/lib/libopenblas.a -lm -fopenmp

OBJS = test_dgemm.o
#OBJS = test_dpotrf.o

all: $(OBJS)

test:
	$(CC) -o test.out $(OBJS) $(LIBS)

run:
	./test.out

clean:
	rm -r *.o
	rm -r test.out
