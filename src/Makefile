CC          := g++

CXXFLAGS        := -O3 -Wall -std=c++11

../bin/hw3: main.cpp
	$(CC) $(CXXFLAGS) -pthread -o $@ $? -fopenmp
clean:
	rm -f ../bin/hw3
