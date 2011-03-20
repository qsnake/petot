atm: libatom
	gfortran -W -Wall -Wextra -fPIC -fbounds-check -fcheck-array-temporaries -o atm.o -c atm.f
	gfortran -o atm atm.o source/libatom.a

libatom:
	cd source; make
