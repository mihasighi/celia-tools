
all: 
	cd src/utils; make all; make install; cd ../../
	cd src/arad/lsum; make all; make install; cd ../../../
	cd src/arad/mset; make all; make install; cd ../../../
	cd src/arad/ucons; make all; make install; cd ../../../
	cd src/shad/sll; make all; make install; cd ../../../
	cd src/shad/shape; make all; make install; cd ../../../

clean:
	cd src/utils; make clean; cd ../../
	cd src/arad/lsum; make clean; cd ../../../
	cd src/arad/mset; make clean; cd ../../../
	cd src/arad/ucons; make clean; cd ../../../
	cd src/shad/sll; make clean; cd ../../../
	cd src/shad/shape; make clean; cd ../../../


tests:
	cd tests/arad/ucons; make all; cd ../../../

