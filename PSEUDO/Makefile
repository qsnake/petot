atm: libatom
	gfortran -o atm.o -c atm.f
	gfortran -o atm atm.o source/libatom.a

libatom:
	cd source; make
