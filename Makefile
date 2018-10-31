all: main_queue

main_queue: main_queue.cc
	g++ -std=c++11 -pthread -O3 main_queue.cc -o main_queue -L/opt/intel/cnc/1.0.100/lib/intel64 -lcnc -lrt -ltbb -ltbbmalloc

clean:
	rm -rf main_queue		
