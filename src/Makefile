CXX = g++
CXXFLAGS = -std=c++17 -Iinclude

SRC = braidgap.cpp lib/braiding.cpp lib/cbraid.cpp
OBJ = $(SRC:.cpp=.o)

TARGET = braidgap
JSONFILE = todo.json

all: $(TARGET) $(JSONFILE)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(JSONFILE):
	touch $(JSONFILE)

clean:
	rm -f $(OBJ) $(TARGET)
	rm -f $(JSONFILE)