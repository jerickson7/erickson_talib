default: libeta.a

libeta.a: cycle_measurement.o filters.o common.o entropy.o volume.o
	ar rcs $@ $^

cycle_measurement.o: cycle_measurement.cpp cycle_measurement.h
	g++ -c -O2 -fPIC $<

filters.o: filters.cpp filters.h
	g++ -c -O2 -fPIC $<

common.o: common.cpp common.h
	g++ -c -O2 -fPIC $<

entropy.o: entropy.cpp entropy.h
	g++ -c -O2 -fPIC $<

volume.o: volume.cpp volume.h
	g++ -c -O2 -fPIC $<

clean:
	rm *.o *.a
