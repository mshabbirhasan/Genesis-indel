CXX = g++ -std=c++11 

genesis:
	$(CXX) src/genesis-indel.cpp -o bin/genesis-indel -lpthread

clean:
	rm -f bin/genesis-indel

.phony: clean default
