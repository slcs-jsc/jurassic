# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = atmfmt brightness climatology day2doy doy2day filter formod hydrostatic interpolate invert jsec2time kernel limb matfmt nadir obs2spec obsfmt planck raytrace retrieval tblfmt tblgen tblrdc time2jsec zenith

# List of tests...
TESTS = cli_test tools_test atm_test mat_test obs_test tbl_test formod_test ret_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib

# Compiler (gcc, icx, nvc, or clang)...
COMPILER ?= gcc

# User-defined preprocessor macros...
DEFINES ?=

# MPI compilation...
MPI ?= 0

# MPI C compiler wrapper...
MPICC ?= mpicc

# Compile for GPUs...
GPU ?= 0

# OpenACC target for GPU builds (gpu, multicore, or host)...
GPU_TARGET ?= gpu

# GPU memory pinning (for GH200)...
GPU_PIN ?= 0

# Static linking...
STATIC ?= 0

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage report...
COV ?= 0

# Path to JURASSIC-UNIFIED library...
UNIDIR = ../../jurassic-unified

# Use JURASSIC-UNIFIED library...
UNIFIED ?= 0

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags 2>/dev/null || echo unknown)

# Set CC and CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -g
ifeq ($(COMPILER),gcc)
  CC = gcc
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
else ifeq ($(COMPILER),icx)
  CC = icx
  CFLAGS += -DHAVE_INLINE -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
else ifeq ($(COMPILER),nvc)
  CC = nvc
  CFLAGS += -mp
else ifeq ($(COMPILER),clang)
  CC = clang
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
else
  $(error Compiler $(COMPILER) not supported)
endif

# Compile for GPUs...
ifeq ($(GPU),1)
  ifeq ($(COMPILER),gcc)
    ifneq ($(GPU_TARGET),gpu)
      $(error GPU_TARGET=$(GPU_TARGET) is not supported by the gcc compiler)
    endif
    CFLAGS += -fopenacc -foffload=-lm
  else ifeq ($(COMPILER),nvc)
    CFLAGS += -acc=$(GPU_TARGET)
    ifeq ($(GPU_TARGET),gpu)
      CFLAGS += -gpu=cc80,cc90
      ifeq ($(GPU_PIN),1)
        CFLAGS += -gpu=mem:separate:pinnedalloc
      endif
    else ifeq ($(GPU_PIN),1)
      $(warning GPU_PIN is ignored for GPU_TARGET=$(GPU_TARGET))
    endif
  else ifeq ($(COMPILER),clang)
    $(error OpenACC is not supported by the clang compiler)
  endif
endif

# LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(STATIC),1)
  LDFLAGS += -lhdf5_hl -lhdf5 -lm -lz -lsz
endif
ifeq ($(GPU),1)
  ifeq ($(COMPILER),nvc)
    ifeq ($(GPU_TARGET),gpu)
      LDFLAGS += -L $(CUDA_PATH)/lib64 -lcudart
    endif
  endif
endif

# Compile with MPI...
ifeq ($(MPI),1)
  CC = $(MPICC)
  CFLAGS += -DMPI
endif

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
  ifeq ($(COMPILER),gcc)
    CFLAGS += -fopt-info
  else ifeq ($(COMPILER),icx)
    CFLAGS += -fopt-info
  else ifeq ($(COMPILER),nvc)
    CFLAGS += -Minfo=all
  else ifeq ($(COMPILER),clang)
    $(error INFO flag not supported by clang compiler)
  endif
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg -fno-omit-frame-pointer
endif

# Compile for coverage...
ifeq ($(COV),1)
  CFLAGS += --coverage
endif

# Static compilation...
ifeq ($(STATIC),1)
  ifeq ($(MPI),1)
    $(error Static compilation does not work for MPI)
  endif
  ifeq ($(UNIFIED),1)
    $(error Static compilation does not work for JURASSIC-UNIFIED)
  endif
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen gpu_test indent install lizard mkdocs mpi_test print-exc strip uninstall $(TESTS)

all: $(EXC)
	rm -f *~

$(EXC): %: %.c jurassic.o $(WRAPPER_OBJ)
	$(CC) $(CFLAGS) -o $@ $< jurassic.o $(WRAPPER_OBJ) $(LDFLAGS)

$(WRAPPER_OBJ): $(WRAPPER).c $(WRAPPER).h Makefile
	$(CC) $(CFLAGS) -c -o  $(WRAPPER_OBJ) $(WRAPPER).c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

check: $(TESTS)

print-exc:
	@printf '%s\n' $(EXC)

$(TESTS) mpi_test gpu_test: all
	@(echo "\n===== Running \"$@\" ... =====") ; \
	  cd ../tests/$@ ; bash ./run.sh \
	  && (echo "\n===== Test \"$@\" passed! =====") \
	  || (echo "\n===== Test \"$@\" failed! =====" ; exit 1)

clean:
	rm -rf $(EXC) *.o *.gcda *.gcno *.gcov coverage* *~

coverage:
	$(MAKE) --no-print-directory clean
	$(MAKE) --no-print-directory COV=1 check
	lcov --capture --directory . --output-file=coverage.info
	genhtml coverage.info --legend --frames --output-directory=coverage

cppcheck:
	cppcheck --enable=all ./

dist:
	cd .. && git ls-files -z | tar --null --files-from=- \
          --transform 's|^|jurassic_$(VERSION)/|' \
          -czvf jurassic_$(VERSION).tgz

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install: all
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
