# C++ compiler (and linker)
CXX=g++
# Compiler flags
CFLAGS=-O3
# Linker flags
LFLAGS=-lpthread -lrt

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=mv

all: $(SRC) $(EXE) $(OBJ)

$(EXE): $(OBJ)
	$(CXX) $(OBJ) $(LFLAGS) -o $@

.cpp.o:
	$(CXX) -c -std=c++0x $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJ) $(EXE)
