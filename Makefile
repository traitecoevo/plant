all:
	make -C doc
	make -C src

clean:
	make -C doc clean
	make -C src clean

test: all
	make -C tests test
