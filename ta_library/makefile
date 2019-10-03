default: libeta

libeta: setup.py etalib.pyx lib/libeta.a
	python3.5 setup.py build_ext --inplace && rm -f etalib.cpp && rm -Rf build

lib/libeta.a:
	make -C lib/

