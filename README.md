transrate-tools
===============

Command-line tools used by [transrate](http://github.com/blahah/transrate) for processing bam files. Written in C++ and using Bamtools.

### building

transrate-tools uses C++11 features, so you'll need at least g++ 4.7 installed. On OSX you can install the latest gcc with `brew install gcc49`. Note that on OSX you always need to tell `cmake` where to find gcc, using the option `-DCMAKE_CXX_COMPILER=$(which g++-4.9)`.

Make sure you clone with submodules:

```bash
$ git clone --recursive git@github.com:cboursnell/transrate-tools.git
```

And you'll need cmake installed.

Next, build bamtools:

on linux:

```bash
$ cd bamtools
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd ../../
```

or on OSX:
```bash
$ cd bamtools
$ mkdir build
$ cd build
$ cmake -DCMAKE_CXX_COMPILER=$(which g++-4.9) ..
$ make
$ cd ../../
```

Then build transrate-tools...

on linux:
```bash
$ cmake .
$ make
```

on OSX:
```bash
$ cmake -DCMAKE_CXX_COMPILER=$(which g++-4.9) .
$ make
```

The executables are called `bam-read` and `bam-split` and will be in the `src` directory.

### bam-split

Parse bam files to separate records that will be filtered by eXpress so that they can be merged back in to the sampled assignments for multi-mapping reads.

### bam-read

Parse bam files after multi-mapping reads have been assigned and aggregate read mapping information by contig.

Currently captured:

 - name
 - bases mapped
 - edit distance
 - bridges
 - length
 - reads mapped
 - both mapped
 - proper pair
 - good
 - uncovered bases
 - mean mapq
 - probability of coverage not being segmented
