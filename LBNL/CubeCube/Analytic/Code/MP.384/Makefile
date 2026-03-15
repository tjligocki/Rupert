CC = gcc

CFLAGS   = -O3 -I./usr/include
CXXFLAGS = -O3
LDFLAGS  = -O3
LIBFLAGS = -L./usr/lib -lgmp -lm

HDRS = matrix.h floatType.h

EXES = iterate1 iterate2 iterate3 bound1 bound2 oneLess root1
OBJS = iterate1.o iterate2.o iterate3.o bound1.o bound2.o \
       oneLess.o root1.o matrix.o

all: $(EXES)

iterate1: iterate1.o floatType.h
	$(CC) $(LDFLAGS) -o iterate1 iterate1.o $(LIBFLAGS)

iterate2: iterate2.o matrix.o matrix.h floatType.h
	$(CC) $(LDFLAGS) -o iterate2 iterate2.o matrix.o $(LIBFLAGS)

iterate3: iterate3.o matrix.o matrix.h floatType.h
	$(CC) $(LDFLAGS) -o iterate3 iterate3.o matrix.o $(LIBFLAGS)

bound1: bound1.o matrix.o matrix.h floatType.h
	$(CC) $(LDFLAGS) -o bound1 bound1.o matrix.o $(LIBFLAGS)

bound2: bound2.o floatType.h
	$(CC) $(LDFLAGS) -o bound2 bound2.o $(LIBFLAGS)

oneLess: oneLess.o floatType.h
	$(CC) $(LDFLAGS) -o oneLess oneLess.o $(LIBFLAGS)

root1: root1.o matrix.o matrix.h floatType.h
	$(CC) $(LDFLAGS) -o root1 root1.o matrix.o $(LIBFLAGS)

$(OBJS): $(HDRS)

clean:
	/bin/rm -f *.o $(EXES)
