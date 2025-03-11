# weaver

*Code is available upon request and will be made publicly available following preprint/publication*

Weaver is a short read mapper for pangenome references. Weaver supports both linear (FASTA) and graphs (rGFA) references inputs and additionally uses small variants in a supplementary VCF file. It is designed to reduce reference bias in variant calling while being a simple drop-in replacement for linear reference mappers.

## Getting started

Go to the [release page](https://github.com/DecodeGenetics/weaver/releases) and get the latest statically linked x86_64 binary. If you prefer, you can also build weaver locally (see below).

### Usage

```sh
# Builds an weaver minimizer index (wmi) at genome.gfa.gz.wmi
weaver index genome.gfa.gz --vcf=small_variants.vcf.gz

# Maps paired short-reads to the graph
weaver map genome.gfa.gz interleaved.fq.gz > out.sam

# Output contains all necessary SAM tags for samtools markdup so it can be piped directly through samtools into a BAM/CRAM
weaver map genome.gfa.gz interleaved.fq.gz | samtools markdup -u - - | samtools view --remove-tag ms -b -o final.bam
samtools index final.bam # Output is already position sorted
```

### Build

Requirements: C++17 compiler (GCC 7+ or clang 12+), CMake 3.2, libzstd, libz.

Recursively git clone the repo and then build with:

```sh
mkdir build && cd build
cmake ..
make weaver
./weaver version # Check version
```

### Test

To run the test suite use:

```sh
make
make test # Runs unit tests
```

## License
MIT License
