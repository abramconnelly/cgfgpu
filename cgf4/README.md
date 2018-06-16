CGF4
===

A tool to manipulate CGF(v4) files.

Compact Genome Format (CGF) files are a compact way to represent whole genomes.
Genomes are represented in a [binary format](../doc/CGFv4-Schema.md) that is a compromise between size
and accessibility.

Whole genomes can be represented in roughly **8.5Mb** to store high quality information.
In order to preserve low quality data, additional space is needed and is dependent on the noise in
the original source.
This can be upwards **30Mb+** but is highly variable.

Quick Start
---

```
$ cgf4 -C mygenome.cgf4
$ cgf4 -e 0 mygenome.cgf4 < test-data/band/hu826751-GS03052-DNA_B01-035e.band 
$ cgf4 -b 0 mygenome.cgf4 
[ 79 8 0 0 0 0 0 -1 0 0 0 389 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 -1 34 -1 185 1]
[ 79 2 0 0 0 0 0 -1 0 0 0 390 0 0 0 0 0 1 0 0 0 0 0 0 26 0 0 1 0 0 -1 34 -1 185 1]
[[ ][ ][ ][ ][ ][ ][ 903 1 ][ ][ 16 1 ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ 96 1 ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ 291 2 ]]
[[ ][ ][ ][ ][ ][ ][ 903 1 ][ ][ 16 1 ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ 96 1 ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ 291 2 ]]
```

Usage
---

```
cgf4 [-H] [-b tilepath] [-e tilepath] [-i ifn] [-o ofn] [-h] [-v] [-V] [ifn]

  [-H|--header]               show header
  [-C|--create-container]     create empty container
  [-I|--info]                 print basic information about CGF file
  [-b|--band tilepath]        output band for tilepath
  [-F fill_level]             bit vector for band fillin (1 canon, 2 cache, 4 ovf, 8 noc, 0xff default)
  [-m|--match]                run concordance on a pair of cgf files (must provide two cgf files)
  [-e|--encode tilepath]      input tilepath band and add it to file, overwriting if it already exists
  [-p|--tilepath tilepath]    tilepath (start)
  [-P|--endtilepath tilepath] end tilepath
  [-s|--tilestep tilestep]    tilestep (start)
  [-S|--endtilestep tilestep] end tilestep
  [--genotype]                genotype flag
  [--all-pairs]               all pairs concordance
  [--show-stats]              show stats
  [-i|--input ifn]            input file (CGF)
  [-o|--output ofn]           output file (CGF)
  [-A|--show-all]             show all tilepaths
  [-h|--help]                 show help (this screen)
  [-v|--version]              show version
  [-V|--verbose]              set verbose level
  [-t|--tilemap tfn]          use tilemap file (instead of default)
  [-Z|--ez-print]             print "ez" structure information
  [-T|--opt-version vopt]     CGF version option.  Must be one of "default" or "noc-inv"
  [-L|--library-version lopt] CGF library version option.  Overwrites default value if specified
  [-U|--update-header]        Update header only
  [-Y|--sanity]               Run sanity checks
  [-q|--hiq]                  Only output high quality information (band output)
```

License
---

AGPLv3

All genomic data used for testing is avaialble under CC0.

Credits
---

Parts of source taken from Curoverse, Inc. ([cgft](https://github.com/curoverse/l7g/tree/master/tools/cgft))


