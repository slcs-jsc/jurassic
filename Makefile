# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = brightness climatology formod hydrostatic interpolate jsec2time kernel limb nadir planck raytrace retrieval time2jsec

# Library directories...
LIBDIR = 

# Include directories...
INCDIR = 

# Linking...
STATIC = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm

# Profiling...
ifdef PROF
  CFLAGS += -O2 -g -pg 
else
  CFLAGS += -O3
endif

# Linking...
ifdef STATIC
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

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

zip:
	zip jurassic_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h
