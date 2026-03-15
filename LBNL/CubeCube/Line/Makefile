all: line sample

line: define.h line.c f.h f.c f1d.h f1d.c brent.h brent.c
	$(CC) -O3 -o line line.c f.c f1d.c brent.c -lm

sample: define.h sample.c f.h f.c
	$(CC) -O3 -o sample sample.c f.c -lm
