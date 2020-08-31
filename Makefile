# Makefile for g++
CC = g++
SRC = main.cpp diagnose.cpp plasma.cpp particles.cpp input.cpp esfield.cpp bfield.cpp 
OBJ := $(SRC:.cpp=.o)
CFLAGS = -fopenmp -pg -O3 -lgsl -lgslcblas -lm -lfftw3 -std=c++11
TARGET = pic

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)
$(OBJ): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm $(OBJ) $(TARGET)
run:$(TARGET)
	-./$(TARGET)

