# weaver

*Source code is available upon request and will be made publicly available prior to publication*

Weaver is a short read mapper for pangenome references. Weaver supports both linear (FASTA) and graphs (rGFA) inputs and additionally uses small variants in a supplementary VCF file. It is designed to reduce reference bias in variant calling while being a simple drop-in replacement for linear reference mappers.

## Getting started

### Installation

Go to the [release page](https://github.com/DecodeGenetics/weaver/releases), download the latest statically linked x86_64 binary and add executable permissions. The binary should run on any 64-bit Linux OS. It does not require a GPU, FPGA or any non-standard hardware. Binary was built and tested on Red Hat Enterprise 9 (RHEL9).

```sh
mkdir -p bin && cd bin
wget https://github.com/DecodeGenetics/weaver/releases/download/v0.2.0/weaver
chmod a+x weaver
```

Optionally you may consider adding the new `bin` directory to your `$PATH`.

### Usage

Weaver has two main subcommands, `weaver index` to construct the Weaver index and `weaver map` to map the reads to the pangenome. The Weaver index is constructed once and loaded during mapping. Typical usage:

```sh
# Builds an weaver minimizer index (wmi) once at genome.gfa.gz.wmi
weaver index genome.gfa.gz --vcf=small_variants.vcf.gz

# Maps paired short-reads to the graph
weaver map genome.gfa.gz interleaved.fq.gz > out.sam

# Output contains all necessary SAM tags for samtools markdup so it can be piped directly through samtools into a BAM/CRAM (requires samtools)
weaver map genome.gfa.gz interleaved.fq.gz | samtools markdup -u - - | samtools view --remove-tag ms -b -o final.bam
samtools index final.bam # Output is already position sorted

# The alignments are in context of the stable sequences of the graph, a linear representation of the GFA. You can make a FASTA file with those sequences with (requires gfatools)
gfatools gfa2fa -s genome.gfa.gz > genome.fa
samtools faidx genome.fa
```

### Build from source code

Requirements: Linux OS (64bit), C++17 compiler (GCC 7+ or clang 12+), CMake 3.2+, libzstd, libz. Tested on RHEL9 OS (64bit).

Recursively git clone the repo and then build with:

```sh
git clone --recursive https://github.com/DecodeGenetics/weaver && cd weaver
mkdir -p build && cd build
cmake .. # Builds a release build by default
make -j4 weaver # Go grab some coffee, building takes a couple of minutes...
```

I recommend to retrying in an empty `build` directory if you encounter any errors from cmake, i.e. some dependency not found.

### Demo

If you want to run a demo with your prebuilt binary or build, you can run the following (working directory should still be your `bin` or `build` directory)

```sh
$ ./weaver version # Check version
0.2.0-5aabbfa
5aabbfac6d86953c78e15e79a36bfcc28d5fa784
$ ./weaver index ../test/data/test_human_10k.gfa.gz --threads 4 -k 17 -w 5 --vcf=../test/data/truth.vcf.gz --log=./demo_weaver_index.log --vverbose
$ ./weaver idxstats ../test/data/test_human_10k.gfa.gz.wmi # Print some index stats
Index stats:
  k = 17
  w = 5
  graph size = 3363 bytes
  map size = 3512 keys
$ ./weaver map ../test/data/test_human_10k.gfa.gz ../test/data/small.read1.fq.gz --fq2=../test/data/small.read2.fq.gz \
    --extra-header-lines=../test/data/extra_header_lines.tsv --threads=4 > ./output.sam
$ grep -v ^@PG ./output.sam | md5sum  # Requires grep and md5sum
7bf29e8b638070e95c9a1313c1eb6f53  -
```

Use the default values for k and w if when running on the full human genome. Details about the available options are in `./weaver [subcommand] --help`

### Unit tests

To run the unit test suite use (possible when source code has been released):

```sh
make # Compiles everything, including the units test
make test # Runs the unit tests
```

## License
MIT License
