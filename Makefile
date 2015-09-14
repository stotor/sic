CXX = g++

#For debugging
#CFLAGS = -Wall -Wextra -Wconversion -g -O0

#For going fast:
CXXFLAGS = -g -Wall -O3

SOURCES = src/main.cpp src/particles.cpp src/utilities.cpp

INCLUDE_FILES = include/particles.hpp include/utilities.hpp

INCLUDE_FOLDER = /Users/stotor/Desktop/1DEM/include

1DEM : $(SOURCES) $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@ -I$(INCLUDE_FOLDER)

clean :
	rm 1DEM
