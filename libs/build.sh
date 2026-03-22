#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in *tar.bz2 ; do tar xvjf $f ; done \
	|| exit

# GSL...
dir=gsl-2.7.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# zlib...
dir=zlib-1.3.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# szip...
dir=szip-2.1.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# HDF5...
dir=hdf5-1.14.4-3
cd $target/src/$dir \
    && ./configure --prefix=$target --with-zlib=$target --with-szlib=$target --enable-hl --disable-fortran --disable-cxx --enable-shared --enable-static \
    && make -j && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.9.2
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --with-plugin-dir=$target/share/netcdf-plugins --disable-dap --disable-byterange --disable-nczarr --disable-libxml2 \
    && make -j && make install && make clean \
	|| exit
