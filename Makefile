all:
	cd src; make

clean:
	cd src; make clean

test: all
	cd tests; make test
