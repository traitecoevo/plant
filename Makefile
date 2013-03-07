all:
	cd src; make

clean:
	cd src; make clean

test: all
	cd tests; Rscript -e "library(testthat); test_dir('.')"
