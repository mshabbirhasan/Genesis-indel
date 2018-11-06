CXX = g++ -std=c++11 

default: all
all: genesis-indel 

genesis:
	$(CXX) $(INCLUDE) src/genesis-indel.cpp -o bin/genesis-indel -lpthread

clean:
	rm -f bin/genesis-indel

.phony: clean default