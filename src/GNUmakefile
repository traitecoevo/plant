## More, and better, warnings
## CXX=c++ # or clang++

## Default GCC for now.
CXX=g++

SOURCES_MODEL  = strategy.cpp parameters.cpp
SOURCES_SPLINE = spline.cpp adaptive_spline.cpp adaptive_spline_r.cpp
SOURCES_ODE    = ode_control.cpp lorenz.cpp ode_r.cpp
SOURCES_R      = interface.cpp
SOURCES_UTIL   = util.cpp

SOURCES = $(SOURCES_MODEL) $(SOURCES_SPLINE) $(SOURCES_ODE) \
	$(SOURCES_UTIL) $(SOURCES_R)

TARGETS = tree.so

all: $(TARGETS)

tree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so
