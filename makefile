CXX=g++
INC_DIR=-I src/
CXXFLAGS=-c -std=c++11 -g -Wall $(INC_DIR)
LDFLAGS=

COMMON_SOURCES= $(wildcard src/*.cc)
COMMON_OBJECTS= $(COMMON_SOURCES:.cc=.o)

TEST_SOURCES= test/test_crystal.cc test/test_distribution.cc
TEST_OBJECTS=   $(TEST_SOURCES:.cc=.o)

TARGET_SOURCES=
TARGET_OBJECTS= $(TARGET_SOURCES:.cc=.o)
TARGET_EXECUTABLE=

TEST_UNITS_SRC= test/units.cc
TEST_UNITS_OBJ= $(TEST_UNITS_SRC:.cc=.o)
TEST_UNITS_EXE= bin/units

COMPARE_STEP_SRC= test/comparisonStep.cc
COMPARE_STEP_OBJ=$(COMPARE_STEP_SRC:.cc=.o)
COMPARE_STEP_EXE= bin/comparisonStep

ALL_OBJ=$(COMMON_OBJECTS) $(TEST_OBJECTS) $(TARGET_OBJECTS) $(TEST_UNITS_OBJ) $(COMPARE_STEP_OBJ)

.PHONY: all clean target units comparisonStep

all: target units comparisonStep

clean:
	rm -f $(ALL_OBJ)
	@echo Clean done.

target: $(TARGET_EXECUTABLE)

units: $(TEST_UNITS_EXE)

comparisonStep: $(COMPARE_STEP_EXE)


$(TARGET_EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TEST_UNITS_EXE): $(COMMON_OBJECTS) $(TEST_OBJECTS) $(TEST_UNITS_OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

$(COMPARE_STEP_EXE): $(COMMON_OBJECTS) $(TEST_OBJECTS) $(COMPARE_STEP_OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

.cc.o:
	$(CXX) $(CXXFLAGS) $< -o $@
