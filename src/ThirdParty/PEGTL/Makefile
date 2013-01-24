# Copyright (c) 2008 Dr. Colin Hirsch
# Please see license.txt for license.

CXX ?= g++

CPPFLAGS := -Iinclude -std=c++0x -pedantic

WARNINGS := -Wall -Wextra
# -Wimplicit -Wconversion -Wcast-align -Woverloaded-virtual -Wold-style-cast -Wformat=2 -Wswitch-enum -Wswitch-default -Wredundant-decls

OPTIMISE := -O1
# -fno-enforce-eh-specs

CXXFLAGS := -march=native $(WARNINGS) $(OPTIMISE)

.PHONY: all clean

SOURCES := $(wildcard */*.cc)
DEPENDS := $(SOURCES:.cc=.d)
PROGRAMS := $(SOURCES:.cc=)

all: $(PROGRAMS)
	@echo "Executing PEGTL unit tests."
	@example/unittests 2> /dev/null

clean:
	rm -f $(DEPENDS) $(PROGRAMS)

.SECONDARY:

%.d: %.cc
	$(CXX) $(CPPFLAGS) -MM -MT $@ $< -o $@

%: %.cc %.d
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

ifeq ($(findstring $(MAKECMDGOALS),clean),)
-include $(DEPENDS)
endif
