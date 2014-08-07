CXX=g++
INC_DIR=-I src/
CXXFLAGS=-c -std=c++11 -g -Wall $(INC_DIR)
LDFLAGS=

COMMON_SOURCES= $(wildcard src/*.cc)
TARGET_SOURCES= compairisonStep.cc
TEST_SOURCES= $(wildcard test/*.cc)

COMMON_OBJECTS= $(COMMON_SOURCES:.cc=.o)
TARGET_OBJECTS= $(TARGET_SOURCES:.cc=.o)
TEST_OBJECTS=   $(TEST_SOURCES:.cc=.o)


TARGET_EXECUTABLE= bin/compairisonStep
TEST_EXECUTABLE= bin/units

.PHONY: all target tests

all: target tests

target: $(TARGET_EXECUTABLE)

tests: $(TEST_EXECUTABLE)

$(TARGET_EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TEST_EXECUTABLE): $(COMMON_OBJECTS) $(TEST_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

.cc.o:
	$(CXX) $(CXXFLAGS) $< -o $@
