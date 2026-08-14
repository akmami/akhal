# AKHAL: Assembly Graph Analysis Tool

## Overview

`akhal` is a toolkit for processing and analyzing r/GFA (Graphical Fragment Assembly) files. 
It validates graphs, reports statistics, and converts GAF (Graph Alignment Format) to SAM (Sequence Alignment Map).

It is structured like `htslib` + `samtools`: a reusable library (`libakhal`) that owns the file formats and data structures, and a thin command-line tool (`akhal`) built on top of it. 

## Installation

Build the library and the CLI with:

```sh
make
```

This produces `libakhal.a` and the `akhal` executable.

Requirements: a C compiler and `zlib` (`-lz`, used by `vg2gfa` to decompress `.vg` files). No protobuf or other external libraries are needed.

## Usage

`akhal` provides the following commands:

```sh
./akhal <PROGRAM> [...ARGS]
```

### Commands

#### 1. `parse`
Validates an r/GFA file and ensures its correctness. It checks segments and links and makes sure that everything is consistent. It also checks the overlapings, if it is presented.

**Usage:**
```sh
./akhal parse <r/GFA file>
```

#### 2. `stats`
Computes and outputs statistics about an r/GFA file.

**Usage:**
```sh
./akhal stats <r/GFA file>
```

The statistics include:
- **Segment count**: Number of segments in the graph.
- **Segment avg length**: Average segment length.
- **Segment std length**: Standard deviation of segment lengths.
- **Segment min length**: Minimum segment length.
- **Segment max length**: Maximum segment length.
- **Link count**: Number of links between segments.
- **Link overlapping avg length**: Average length of overlapping links.
- **Link overlapping std length**: Standard deviation of link overlap lengths.
- **Minimum in degree**: Minimum number of incoming links.
- **Maximum in degree**: Maximum number of incoming links.
- **Minimum out degree**: Minimum number of outgoing links.
- **Maximum out degree**: Maximum number of outgoing links.

#### 3. `extract`
Extract information from the r/GFA file.

**Usage:**
```sh
./akhal gaf2sam extract [OPTION] <r/GFA file> <OUTPUT file>
```

Options:
- **fa**: Reference genome. Output file should end with `.fa` or `.fasta`

#### 4. `sort`
Topologically sorts a graph and renumbers its nodes `1..N` in the sorted order. Ordering uses Kahn's algorithm; ties in the ready set (the "hops") are broken alphabetically by node **sequence content**, so the result is independent of the input's node numbering — two graphs identical in topology and content sort the same way. Nodes in cycles, if any, are appended after the acyclic prefix. The sorted graph is re-emitted with the new ids (S/L/P lines remapped, link orientation and overlap preserved).

**Usage:**
```sh
./akhal sort <input .gfa file> [output .gfa file]
```

Note: If no output file is given, the sorted GFA is written to standard output.

#### 5. `vg2gfa`
Converts vg's native `.vg` format (gzip/BGZF-compressed Protobuf) to GFA. The `.vg` parser is hand-written pure C - no protobuf or libvgio needed, only zlib.

**Usage:**
```sh
./akhal vg2gfa <input .vg file> [output .gfa file]
```

Note: If no output file is given, the GFA is written to standard output.

#### 6. `gaf2sam`
Converts a GAF file to a SAM file.

**Usage:**
```sh
./akhal gaf2sam <r/GFA file> <GAF file>  <FASTA file> <OUTPUT file>  [--simple]
```

Note: The reads should be stored in FASTA format and provided to the program. 
GAF file does not store sequences, hence, reads are needed when converting to SAM.

Note: `simple` option is optional. 
If provided, CIGAR string matches `(=)` and mismatches `(X)` will be replaced with sequence match `(M)`.

#### 7. `sampoke`
Validate SAM file (converted from gaf). It takes reference file and SAM to process CIGAR strings. 
Optionally, it can print the filtered SAM file that contains valid lines.

**Usage:**
```sh
./akhal sampoke <FASTA file> <GAF file> <OUTPUT file>
```

Note: Output file here is optional.

#### 8. `annotate`
Traces the origin of every graph node and saves the result as a binary `.annot` file. 
The backbone reference path (a `P` line) is identified first - by `--ref <name>`, or defaulting to the graph's first path - and its nodes are marked `backbone` with their reference coordinates. 
With `--vcf`, each variant is matched to the alternative side of its bubble: the shared REF/ALT prefix is stripped, the branch point is located on the backbone, and the nodes spelling the alternate allele are annotated like `SNP chr1:12345 A>G rs99`. 
Pure deletions produce only an edge, so they have nothing to annotate. 
With `--fasta`, each sequence is walked through the graph (link overlaps honoured) and every visited non-backbone node is annotated `SEQ <name> <offset>` with its origin and position in that sequence. 
Both inputs are optional; nodes that no input explains are simply left out and later report as `unknown`. 
A node explained more than once accumulates its annotations separated by `; `.

**Usage:**
```sh
./akhal annotate <r/GFA file> <OUTPUT .annot file> [--vcf <VCF file>] [--fasta <FASTA file>] [--ref <path name>]
```

#### 9. `annotget`
Looks up node annotations in a `.annot` file - the graph itself is not needed. 
Each queried node prints one line, `<id> <status> [<info>]`, where status is `backbone`, `annot`, or `unknown`. 
With no node ids, every record in the file is dumped.

**Usage:**
```sh
./akhal annotget <.annot file> [node id ...]
```

Example:
```sh
$ ./akhal annotget graph.annot 4 6 2 999
4	annot	SNP chr1:5 A>C rs1; SEQ sample1 4
6	annot	INS chr1:9 C>CTT
2	backbone	REF chr1 4-5
999	unknown
```

## License
is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. 
For more detailed terms, please refer to the [license file](https://github.com/akmami/akhal/blob/main/LICENCE).

## Author
Developed by Akmuhammet Ashyralyyev.

**Note from the author**:

This tool is named after one of the most elegant horses, the [Akhal-Teke](https://en.wikipedia.org/wiki/Akhal-Teke). 
This breed is one of the oldest domesticated animals and is considered one of the most beautiful and intelligent horses in the world.
