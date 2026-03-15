CFLAGS = -O3 -I${HOME}/Software/Linux/include
CXXFLAGS = -O3 -I${HOME}/Software/Linux/include
LDFLAGS = -O3

LIBS = -L${HOME}/Software/Linux/lib -lmpfr -lgmp -lm

HDRS = rupert.h

EXES = check1 check2 improve1 gradient \
       rotateRelax study1.002.003 gradientRelax

OBJS = check1.o check2.o improve1.o gradient.o \
       rotateRelax.o study1.002.003.o gradientRelax.o

LIBOBJS = rupert.o

all: $(EXES)

check1: check1.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o check1 check1.o $(LIBOBJS) $(LIBS)

check2: check2.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o check2 check2.o $(LIBOBJS) $(LIBS)

improve1: improve1.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o improve1 improve1.o $(LIBOBJS) $(LIBS)

gradient: gradient.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o gradient gradient.o $(LIBOBJS) $(LIBS)

gradientRelax: gradientRelax.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o gradientRelax gradientRelax.o $(LIBOBJS) $(LIBS)

rotateRelax: rotateRelax.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o rotateRelax rotateRelax.o $(LIBOBJS) $(LIBS)

study1.002.003: study1.002.003.o $(LIBOBJS)
	$(CXX) $(LDFLAGS) -o study1.002.003 study1.002.003.o $(LIBOBJS) $(LIBS)

$(OBJS): $(HDRS)

$(LIBOBJS): $(HDRS)
