#! /bin/sh
g++ -Wl,-Bstatic -std=c++11 transrate-bam-split.cpp -o bam-split -L bamtools/lib -lbamtools -lbamtools-utils -static-libgcc -I bamtools/include -lz
