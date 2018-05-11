all: main

main: main.cc
	g++ -std=c++11 -pthread -O3 main.cc -o main -L/opt/intel/cnc/1.0.100/lib/intel64 -lcnc -lrt -ltbb -ltbbmalloc

clean:
	rm -rf main		
