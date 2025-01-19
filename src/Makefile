# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = brightness climatology filter formod hydrostatic interpolate invert jsec2time kernel limb nadir obs2spec planck raytrace retrieval tblfmt tblgen time2jsec

# List of tests...
TESTS = limb_test nadir_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib

# Path to JURASSIC-UNIFIED library...
UNIDIR = ../../jurassic-unified

# Use JURASSIC-UNIFIED library...
UNIFIED ?= 0

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage report...
COV ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)

# Compiler...
CC ?= gcc

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp $(OPT) -g

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm

# Use JURASSIC-UNIFIED library...
ifeq ($(UNIFIED),1)

  # Wrapper file...
  WRAPPER=$(UNIDIR)/include/jurassic_unified_library
  WRAPPER_OBJ=wrapper.o

  # Add paths...
  INCDIR += -I $(UNIDIR)/include -I .
  LIBDIR += -L $(UNIDIR)/unified_library -L $(CUDA_PATH)/lib64

  # Add flags...
  CFLAGS += -DUNIFIED
  LDFLAGS += -lstdc++ -lcudart -ljurassic_unified

endif

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Compile for coverage...
ifeq ($(COV),1)
  CFLAGS += --coverage
endif

# Static compilation...
ifeq ($(STATIC),1)
  ifeq ($(UNIFIED),1)
    $(error Static compilation does not work for JURASSIC-UNIFIED)
  else
    CFLAGS += -static
  endif
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen indent \
           install lizard mkdocs strip uninstall

all: $(EXC)
	rm -f *~

$(EXC): %: %.c jurassic.o $(WRAPPER_OBJ)
	$(CC) $(CFLAGS) -o $@ $< jurassic.o $(WRAPPER_OBJ) $(LDFLAGS)

$(WRAPPER_OBJ): $(WRAPPER).c $(WRAPPER).h Makefile
	$(CC) $(CFLAGS) -c -o  $(WRAPPER_OBJ) $(WRAPPER).c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

check: $(TESTS)

$(TESTS): all
	@(echo "\n===== Running \"$@\" ... =====") ; \
	  cd ../tests/$@ ; ./run.sh \
	  && (echo "\n===== Test \"$@\" passed! =====") \
	  || (echo "\n===== Test \"$@\" failed! =====" ; exit 1)

clean:
	rm -rf $(EXC) *.o *.gcda *.gcno *.gcov coverage* *~

coverage:
	lcov --capture --directory . --output-file=coverage.info ; \
	genhtml coverage.info --legend --frames --output-directory=coverage ; \
	gcov -abcfu *.c ; echo "FINISHED" || true

cppcheck:
	cppcheck --enable=all ./

dist:
	cd ..  && tar cvfz jurassic_$(VERSION).tgz bin/* \
          docs/Doxyfile libs/build.sh libs/*bz2 \
          projects/limb/* projects/nadir/* src/* tests/*

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
