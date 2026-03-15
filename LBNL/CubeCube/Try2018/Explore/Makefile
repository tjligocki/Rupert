CC = gcc
CFLAGS = -O3

all: makeIdentity minimize explore explore1 explore2 explore3 exploreRandom \
     exploreLine exploreLines exploreMinimum minimizeLines exploreRotation  \
     rowGrads planeRotation

tools.o: tools.h

makeIdentity: makeIdentity.o tools.o
	$(CC) $(CFLAGS) -o makeIdentity makeIdentity.o tools.o -lm

minimize: minimize.o tools.o
	$(CC) $(CFLAGS) -o minimize minimize.o tools.o -lm

explore: explore.o tools.o
	$(CC) $(CFLAGS) -o explore explore.o tools.o -lm

explore1: explore1.o tools.o
	$(CC) $(CFLAGS) -o explore1 explore1.o tools.o -lm

explore2: explore2.o tools.o
	$(CC) $(CFLAGS) -o explore2 explore2.o tools.o -lm

explore3: explore3.o tools.o
	$(CC) $(CFLAGS) -o explore3 explore3.o tools.o -lm

exploreRandom: exploreRandom.o tools.o
	$(CC) $(CFLAGS) -o exploreRandom exploreRandom.o tools.o -lm

exploreLine: exploreLine.o tools.o
	$(CC) $(CFLAGS) -o exploreLine exploreLine.o tools.o -lm

exploreLines: exploreLines.o tools.o
	$(CC) $(CFLAGS) -o exploreLines exploreLines.o tools.o -lm

exploreMinimum: exploreMinimum.o tools.o
	$(CC) $(CFLAGS) -o exploreMinimum exploreMinimum.o tools.o -lm

minimizeLines: minimizeLines.o tools.o
	$(CC) $(CFLAGS) -o minimizeLines minimizeLines.o tools.o -lm

exploreRotation: exploreRotation.o tools.o
	$(CC) $(CFLAGS) -o exploreRotation exploreRotation.o tools.o -lm

rowGrads: rowGrads.o tools.o
	$(CC) $(CFLAGS) -o rowGrads rowGrads.o tools.o -lm

planeRotation: planeRotation.o tools.o
	$(CC) $(CFLAGS) -o planeRotation planeRotation.o tools.o -lm
