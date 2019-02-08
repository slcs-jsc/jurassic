#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)
threads=$(cat /proc/cpuinfo | grep processor | wc -l)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
    || exit

# GSL...
dir=gsl-2.5
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$threads && make check && make install && make clean \
    || exit
