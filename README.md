transrate-bam-read
==================

Command-line tool using the bamtools c++ library to parse bam files and aggregate read mapping
information by contig for [transrate](https://github.com/Blahah/transrate).

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

## Building

Make sure you clone with submodules:

```bash
$ git clone --recursive git@github.com:cboursnell/transrate-bam-read.git
```

Then just run the build script:

```bash
$ ./build
```

The executable is called `bam-read`.