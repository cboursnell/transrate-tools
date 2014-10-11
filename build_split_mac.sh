#! /bin/sh

g++-4.9 -std=c++11 -o bam-split bamtools/lib/libbamtools.a -lbamtools -I bamtools/include -lz transrate-bam-split.cpp
