# Makefile for Molecular Trajectory Analysis

# Compiler settings - Using NVIDIA HPC SDK
CXX = nvc++
CXXFLAGS = -O3 -std=c++17 -Wall

# Target executable
TARGET = trajectory_analyzer

# Source files
SRCS = main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = config.h trajectory.h cluster.h

# Default target
all: $(TARGET)

# Link the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Run the program with a sample config file
run: $(TARGET)
	./$(TARGET) config.txt

.PHONY: all clean run
