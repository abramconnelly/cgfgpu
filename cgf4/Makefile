CC=g++
#CC_FLAGS=-O3 -std=c++11 -msse4.2
CC_FLAGS=-g -std=c++11 -msse4.2
CC_LIB=-lsdsl -L$(HOME)/lib
CC_INC=-I$(HOME)/include

src=$(wildcard *.cpp)
obj=$(src:.cpp=.o)

%.o: %.cpp
	$(CC) -c -o $@ $^ $(CC_FLAGS) $(CC_LIB) $(CC_INC)

gcgf: $(obj)
	$(CC) $(CC_FLAGS) $^ -o cgf4 $(CC_LIB) $(CC_INC)

.PHONY: clean
clean:
	rm -f *.o
	rm -f cgf4

