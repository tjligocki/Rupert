# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += equalRows

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = equalRows

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../../../../Chombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = 

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
