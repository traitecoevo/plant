SOURCES_MODEL  = strategy.cpp control.cpp parameters.cpp \
	disturbance.cpp environment.cpp \
	plant.cpp cohort_discrete.cpp plant_spline.cpp plant_approx.cpp \
	species.cpp patch.cpp cohort_top.cpp \
	cohort_schedule.cpp ebt.cpp metacommunity.cpp
SOURCES_SPLINE = spline.cpp adaptive_spline.cpp adaptive_spline_r.cpp \
	multi_spline.cpp
SOURCES_INTEGRATION = quadrature.cpp quadrature_rules.cpp \
	integration.cpp integration_workspace.cpp
SOURCES_ODE    = ode_target.cpp ode_control.cpp lorenz.cpp ode_r.cpp
SOURCES_R      = interface.cpp
SOURCES_UTIL   = util.cpp lookup.cpp find_root.cpp functor.cpp \
	integrator.cpp gradient.cpp

SOURCES = $(SOURCES_MODEL) $(SOURCES_SPLINE) $(SOURCES_INTEGRATION) \
	$(SOURCES_ODE) $(SOURCES_UTIL) $(SOURCES_R)

TARGETS = tree.so

all: $(TARGETS)

tree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so
	rm -rf tree.so.dSYM

.PHONY: all clean
