transrate-tools
===============

Command-line tools used by [transrate](http://github.com/blahah/transrate) for processing bam files. Written in C++ and using Bamtools.

### building

Make sure you clone with submodules:

```bash
$ git clone --recursive git@github.com:cboursnell/transrate-tools.git
```

And you'll need cmake installed.

Then just...

```bash
$ cmake .
$ make
```

On OSX, you will need to have gcc-4.7 or later installed (you can do this with `brew install gcc49`), and will need to pass the location of gcc to cmake:

```bash
$ cmake -DCMAKE_CXX_COMPILER=$(which g++-4.9) .
$ make
```

The executables are called `bam-read` and `bam-split`.

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
