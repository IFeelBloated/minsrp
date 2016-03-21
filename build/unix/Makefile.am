warningflags = -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers -Wshadow
includeflags = -I$(srcdir)/../../src
commonflags = -O3 $(MFLAGS) $(warningflags) $(includeflags)
AM_CXXFLAGS = -std=c++11 $(commonflags)

lib_LTLIBRARIES = libminsrp.la

libminsrp_la_SOURCES = ../../src/Kernel.hpp \
                        ../../src/Source.cpp

libminsrp_la_LDFLAGS = -no-undefined -avoid-version $(PLUGINLDFLAGS)
