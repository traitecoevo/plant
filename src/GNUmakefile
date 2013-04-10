SOURCES_MODEL  = strategy.cpp parameters.cpp plant.cpp \
	cohort.cpp cohort_discrete.cpp patch.cpp
SOURCES_SPLINE = spline.cpp adaptive_spline.cpp adaptive_spline_r.cpp
SOURCES_ODE    = ode_target.cpp ode_control.cpp lorenz.cpp ode_r.cpp
SOURCES_R      = interface.cpp
SOURCES_UTIL   = util.cpp lookup.cpp find_root.cpp functor.cpp \
	integrator.cpp

SOURCES = $(SOURCES_MODEL) $(SOURCES_SPLINE) $(SOURCES_ODE) \
	$(SOURCES_UTIL) $(SOURCES_R)

TARGETS = tree.so

# TODO: To work around issues of dependencies, we can use 
#    R CMD COMPILE $(FILE)
# as the .cpp->.o rule, and I can then have files depend on each
# other, and on other files.
#
# However, this still uses an internal rule set for determining if the
# file needs recompilation, so we'd need a rm -f $(FILE).o first.

all: $(TARGETS)

tree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so
	rm -rf tree.so.dSYM
