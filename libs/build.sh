#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in *.tar.bz2 ; do tar xvjf $f ; done \
    || exit

# GSL...
dir=gsl-2.7.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
    || exit
