## More, and better, warnings
## CXX=c++ # or clang++

## Default GCC for now.
CXX=g++

SOURCES_SPLINE = Spline.cpp AdaptiveSpline.cpp RAdaptiveSpline.cpp
SOURCES_ODE    = Step.cpp Control.cpp Evolve.cpp Lorenz.cpp
SOURCES_R      = interface.cpp
SOURCES_UTIL   = util.cpp

SOURCES = $(SOURCES_SPLINE) $(SOURCES_ODE) $(SOURCES_UTIL) $(SOURCES_R)

TARGETS = tree.so

all: $(TARGETS)

tree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so
