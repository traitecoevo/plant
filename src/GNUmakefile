## More, and better, warnings
## CXX=c++ # or clang++

## Default GCC for now.
CXX=g++

SOURCES = util.cpp          \
	Spline.cpp          \
	AdaptiveSpline.cpp  \
	RAdaptiveSpline.cpp \
	interface.cpp

TARGETS = tree.so

all: $(TARGETS)

tree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so
