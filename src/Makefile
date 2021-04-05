# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = brightness climatology formod hydrostatic interpolate invert jsec2time kernel limb nadir planck raytrace retrieval time2jsec

# Library directories...
LIBDIR ?= -L ../libs/build/lib

# Include directories...
INCDIR ?= -I ../libs/build/include

# Static compilation...
STATIC ?= 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
CC ?= gcc

# CFLAGS...
CFLAGS ?= $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -O3 -g

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm

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

clean:
	rm -f $(EXC) *.o *~

dist:
	zip jurassic_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	cd ../docs && doxygen && cd ../docs/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC)
