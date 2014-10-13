#! /bin/sh
mkdir bamtools/build
cd bamtools/build
cmake ../
make
cd ../..
g++ -std=c++11 segmenter.cpp transrate_pileup.cpp transrate-bam-read.cpp -o bam-read -lbamtools bamtools/lib/libbamtools.a -I bamtools/include -lz
