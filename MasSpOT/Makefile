optimal_transport.so: optimal_transport.o
	g++ -shared -g -o optimal_transport.so optimal_transport.o

optimal_transport.o: optimal_transport.cpp optimal_transport.h
	g++ -c -Wall -g -O3 -msse2  -fpic -I eigen optimal_transport.cpp
