#CXX = g++
CXX = mpicxx

#For debugging
#CFLAGS = -Wall -Wextra -Wconversion -g -O0

#For going fast:
CXXFLAGS = -Wall -O3

SOURCES = src/main.cpp src/particlespecies.cpp \
	  src/speciesgroup.cpp src/fields.cpp src/utilities.cpp 

INCLUDE_FILES = include/particlespecies.hpp include/fields.hpp \
	        include/speciesgroup.hpp include/utilities.hpp

INCLUDE_FOLDER = include

sic : $(SOURCES) $(INCLUDE_FILES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@ -I$(INCLUDE_FOLDER)

clean :
	rm sic
