
all: 
	cd src/utils; make all; make install; cd ../../
	cd src/arad/ucons; make all; make install; cd ../../../

tests:
	cd tests/arad/ucons; make all; cd ../../../

