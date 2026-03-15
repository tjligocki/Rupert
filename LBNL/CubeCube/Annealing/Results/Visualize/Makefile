CFLAGS   = -O3
CXXFLAGS = -O3
LDFLAGS  = -O3

EXES = case12 case13 case23 neighbor23 neighbor34 neighbor neighbor2

all: $(EXES)

case12: case12.o
	$(CXX) $(LDFLAGS) -o case12 case12.o -lm

case13: case13.o
	$(CXX) $(LDFLAGS) -o case13 case13.o -lm

case23: case23.o
	$(CXX) $(LDFLAGS) -o case23 case23.o -lm

neighbor23: neighbor23.o
	$(CXX) $(LDFLAGS) -o neighbor23 neighbor23.o -lm

neighbor34: neighbor34.o
	$(CXX) $(LDFLAGS) -o neighbor34 neighbor34.o -lm

neighbor: neighbor.o
	$(CXX) $(LDFLAGS) -o neighbor neighbor.o -lm

neighbor2: neighbor2.o
	$(CXX) $(LDFLAGS) -o neighbor2 neighbor2.o -lm

clean:
	/bin/rm -f *.o $(EXES)
