all:
	cd src; make -j9
	cp -v src/*.x .

clean:
	cd src; make clean
	rm -v *.x