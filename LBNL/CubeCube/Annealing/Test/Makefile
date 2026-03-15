MATH = -lm

LIBS = $(MATH)

CFLAGS = -O3
CC = gcc

LDFLAGS	= -O3

HFILES = define.h global.h
CFILES = cube.c params.c new.c copy.c initchain.c next.c energy.c	\
	 approve.c printchain.c debugchain.c panic.c check.c orthonormalize.c
OFILES = cube.o params.o new.o copy.o initchain.o next.o energy.o	\
	 approve.o printchain.o debugchain.o panic.o check.o orthonormalize.o

all: cube

cube: $(OFILES)
	$(CC) $(LDFLAGS) -o cube $(OFILES) $(LIBS)

$(OFILES): $(HFILES)

clean:
	/bin/rm -f *.o cube
