# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = brightness climatology filter formod hydrostatic interpolate invert jsec2time kernel limb nadir obs2spec planck raytrace retrieval tblfmt tblgen time2jsec

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Compiler...
CC ?= gcc

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp $(OPT) -g

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c jurassic.o
	$(CC) $(CFLAGS) -o $@ $< jurassic.o $(LDFLAGS)

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

check: all
	cd ../example/nadir && ./run.sh
	cd ../example/limb && ./run.sh

clean:
	rm -f $(EXC) *.o *~

dist:
	zip jurassic_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	cd ../docs && doxygen && cd ../docs/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install: all
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

strip:
	strip $(EXC)
